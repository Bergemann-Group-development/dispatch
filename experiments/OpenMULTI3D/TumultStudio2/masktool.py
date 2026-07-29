"""The LineMasks tool: import, view and edit the shared set of fit/synth windows.

Toolbar **Masks** opens a right-sidebar panel that imports a linemask file (``input/linemasks``),
lists the masks, tunes the synth-window padding, toggles per-mask local renormalisation, and
saves the edited set back out. The masks themselves live in the shared :class:`maskstore.MaskStore`
(``ctx.masks``) so the canvas — and, later, the abundance / synthesize-fit workflows — all read one
source of truth; every mutation broadcasts ``"masks-changed"`` and this tool repaints from it.

The canvas draws the set while the tool is active: each mask's synth/visual window as a faint red
band, its fit window amber, and the active mask accented. Editing the windows by dragging (box
select) and the local-continuum markers are layered on in later phases.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

import astropy.units as u
import numpy as np
from bokeh.models import TextInput

from canvas import ACCENT
from controls import Field, Form
from filepicker import FilePicker
from htmlbridge import parse_action, pulse
from linemask import seed_cont_markers, write_linemask_file
from rangeselect import RangeSelect
from tool import Tool
from units import to_range


class MaskTool(Tool):
    name = "linemask"

    def __init__(self, ctx, store, input_dir):
        super().__init__(ctx)
        self.store = store                      # shared MaskStore (also ctx.masks)
        self.input_dir = str(input_dir)
        self._renorm_uid = None                 # spectrum currently carrying the local renorm
        self.own_panel("linemask-panel-content")   # sidebar: import/pad/save controls
        self._table_bridge = ctx.bridges.html("linemask-table")   # centre: the mask table
        # the sidebar panel, the toolbar button, AND the centre bottom-panel (mask table below
        # the spectrum, like the abundance tool's #abund-bottom) all toggle together.
        self.own_mode_bridges((".sidebar-right", "linemask-active"),
                              ("#linemask-btn", "active"),
                              ("#main", "linemask-mode"))

        # synth-window padding (Å): re-derives every mask's [smin,smax] from its fit window
        self.pad_form = Form("linemask", [
            Field("pad", "Synth padding [Å]", kind="number", value=store.pad, min=0.0, step=0.1),
        ], on_change=self._on_pad)

        # inline import picker, anchored at input/linemasks (reveal toggled from the panel)
        self.open_bridge = ctx.bridges.cls("#linemask-picker", "open")
        self.picker = FilePicker(self._on_import, container_id="linemask-results",
                                 root=input_dir, start_dir=input_dir)
        self.picker.layout.name = "linemask_picker"

        # per-row / save / clear actions come back on a dedicated space-delimited channel
        self.actions = TextInput(value="", visible=False, name="linemask-actions")
        self.actions.on_change("value", self._on_action)

        # box-select gesture to place / resize a mask window (rectangular select, per the user's
        # preference over edge-drag handles); armed only while the tool owns the sidebar
        self.range = RangeSelect(self.canvas, ACCENT, "Draw a line-mask window")
        self.range.on_select = self._on_range

        ctx.subscribe("masks-changed", lambda **_: self._on_masks_changed())

    # --- doc wiring -----------------------------------------------------------------------
    @property
    def models(self):
        return [self.actions, *self.pad_form.models]

    # --- Tool protocol --------------------------------------------------------------------
    def activate(self):
        self.active = True
        self._render_canvas()
        self._apply_renorm()
        self._render_panels()
        self.range.arm()             # box-select places / resizes mask windows
        self.push_mode("on")

    def deactivate(self):
        if not self.active:
            return
        self.active = False
        self.range.disarm()
        self.range.clear_selection()
        self.canvas.clear_masks()
        self.push_mode("off")

    def update(self):
        """Canvas ``"changed"`` (unit flip, selection change …): repaint the bands and re-apply
        the local renorm (which follows the selected spectrum)."""
        if self.active:
            self._render_canvas()
            self._apply_renorm()

    # --- reactions ------------------------------------------------------------------------
    def _on_masks_changed(self):
        self.pad_form.fields["pad"].value = self.store.pad
        self._render_panels()
        if self.active:
            self._render_canvas()
            self._apply_renorm()

    def _render_canvas(self):
        self.canvas.render_masks(self.store.masks, self.store.active)

    def _render_panels(self):
        """Repaint the sidebar controls and the centre mask table."""
        self.push_panel()
        self._table_bridge.value = self.table_html()

    def _apply_renorm(self):
        """Fold every renorm-on mask's local linear continuum into the *selected* spectrum's
        continuum (injected spline knots) and redraw it. The renorm follows the selection, so if
        it moved, the previous target is reverted first."""
        obj = self.canvas.selected()
        # seed markers (a percentile linear fit) for any renorm mask that lacks them
        if obj is not None:
            for m in self.store.masks:
                if m.renorm and not m.cont_pts:
                    m.cont_pts = seed_cont_markers(obj.raw, m)
        windows = [(m.smin, m.smax, m.renorm_knots()) for m in self.store.masks
                   if m.renorm and m.renorm_knots()]
        prev = self._renorm_uid
        if prev and (obj is None or obj.uid != prev):
            old = self.canvas._find(prev)
            if old is not None:
                old.renorm_windows = []
                self.canvas.renormalise(old)
            self._renorm_uid = None
        if obj is not None:
            obj.renorm_windows = windows
            self.canvas.renormalise(obj)
            self._renorm_uid = obj.uid if windows else None

    def _jump_to_active(self):
        """Zoom the spectrum to the active mask's window so clicking a table row reveals it."""
        m = self.store.active_mask()
        if m is None:
            return
        lo, hi = (m.lmin, m.lmax)
        if not (np.isfinite(lo) and np.isfinite(hi)):
            lo, hi = m.smin, m.smax          # fall back to the synth window (e.g. no fit bounds)
        self.canvas.focus_range(lo, hi)

    def _on_pad(self, key, value):
        try:
            self.store.set_pad(float(value))
        except (TypeError, ValueError):
            pass

    def _on_range(self, x0, x1):
        """A box was drawn: resize the active mask's fit window, or create a new mask when none
        is active. ``x0``/``x1`` are in the canvas display unit; masks store Å."""
        lo, hi = to_range(x0, x1, self.canvas.xunit, u.AA)
        if self.store.active_mask() is not None:
            self.store.set_fit(self.store.active, lo, hi)
            self.notify(f"Resized mask to {lo:.3f}–{hi:.3f} Å.")
        else:
            self.store.add(lo, hi, lam0=0.5 * (lo + hi))
            self.notify(f"Added mask {lo:.3f}–{hi:.3f} Å.")
        self.range.clear_selection()   # the mask band now shows the window

    def _on_import(self, path):
        try:
            mset = self.store.load_file(path)
        except Exception as exc:
            logging.exception("Line-mask import failed")
            self.notify(f"Could not read line masks: {exc}")
            return
        pulse(self.open_bridge, "off")     # collapse the inline picker
        self.notify(f"Imported {len(mset)} line masks from {Path(path).name}.")

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "select":
            self.store.set_active(_int(arg))
            self._jump_to_active()
        elif kind == "del":
            self.store.remove(_int(arg))
        elif kind == "renorm":
            i = _int(arg)
            m = self.store.masks[i] if 0 <= i < len(self.store.masks) else None
            if m is not None:
                self.store.set_renorm(i, not m.renorm)
        elif kind == "new":
            self.store.set_active(None)    # next box-select creates a fresh mask
            self.notify("Draw a box on the spectrum to place the new mask.")
        elif kind == "clear":
            self.store.clear()
        elif kind == "save":
            self._save()

    def _save(self):
        if not self.store.masks:
            self.notify("No line masks to save — import a file first.")
            return
        src = self.store.mset.source_path
        out = (Path(src).with_name(Path(src).stem + "_edited.txt") if src
               else Path.cwd() / "linemask_export.txt")
        try:
            write_linemask_file(out, self.store.mset)
        except Exception as exc:
            logging.exception("Line-mask save failed")
            self.notify(f"Could not save line masks: {exc}")
            return
        self.notify(f"Saved {len(self.store.masks)} line masks to {out.name}.")

    # --- rendering ------------------------------------------------------------------------
    def panel_html(self):
        """Sidebar controls: synth padding, a count/source line, and save/clear (the mask table
        itself lives in the centre bottom panel, see :meth:`table_html`)."""
        n = len(self.store.masks)
        src = self.store.mset.source_path
        srcline = (f'<div class="hint">Loaded {n} mask{"" if n == 1 else "s"} from '
                   f'{html.escape(Path(src).name)}.</div>' if src and n else
                   '<div class="hint">Import a line-mask file to load fit windows.</div>')
        active = self.store.active_mask() is not None
        new_hint = ("Draw a box to resize the active mask." if active else
                    "Draw a box on the spectrum to create a mask.")
        new_btn = (f'<button class="text-btn lm-new" title="{new_hint}">+ New mask</button>')
        saveclear = ('<div class="lm-actions">'
                     '<button class="text-btn lm-save">Save masks…</button>'
                     '<button class="text-btn lm-clear">Clear</button></div>') if n else ""
        return (self.pad_form.html()
                + f'<div class="section-header">{n} line mask{"" if n == 1 else "s"}</div>'
                + srcline + new_btn + f'<div class="hint">{new_hint}</div>' + saveclear)

    def table_html(self):
        """The centre bottom-panel mask table: one row per mask, click to make it active."""
        masks = self.store.masks
        if not masks:
            return ('<div class="lm-empty">No line masks loaded — use '
                    '<b>Import masks…</b> in the Masks panel.</div>')
        rows = []
        for i, m in enumerate(masks):
            sel = " lm-active" if i == self.store.active else ""
            lam0 = f"{m.lam0:.3f}" if m.lam0 == m.lam0 else "—"   # NaN check
            rn = " checked" if m.renorm else ""
            rows.append(
                f'<tr class="lm-row{sel}" data-idx="{i}">'
                f'<td>{lam0}</td>'
                f'<td>{m.lmin:.3f} – {m.lmax:.3f}</td>'
                f'<td>{m.smin:.3f} – {m.smax:.3f}</td>'
                f'<td class="lm-note">{html.escape(m.note)}</td>'
                f'<td><input type="checkbox" class="lm-renorm" data-idx="{i}"{rn} '
                f'title="Local renormalisation"></td>'
                f'<td><button class="icon-btn lm-del" data-idx="{i}" title="Remove">'
                f'<i class="ti ti-x"></i></button></td></tr>')
        return (
            '<table class="linelist-table lm-table">'
            '<thead><tr><th>λ₀ [Å]</th><th>fit window [Å]</th><th>synth window [Å]</th>'
            '<th>note</th><th title="Local renormalisation">rn</th><th></th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


def _int(arg, default=-1):
    try:
        return int(arg)
    except (TypeError, ValueError):
        return default
