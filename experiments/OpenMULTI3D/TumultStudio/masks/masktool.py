"""The LineMasks view: import, view and edit the sets of fit/synth windows.

A line-mask set comes from the toolbar's **Open → Line mask set…** (read from a file) or
**New → Line mask set…** (started empty) commands, one pane per :class:`masksources.MaskSource`
in the same modal; generation from a reference spectrum is a third, currently disabled source. However it was born it becomes a *line-mask object* listed in
the left sidebar under its name, carrying the element it diagnoses. Clicking one of those objects opens this view: the right sidebar gets the
mask controls (synth-window padding, save-as / revert) and the centre gets the mask table below
the spectrum, with every window drawn on the canvas.

The sets themselves live in the shared :class:`maskstore.MaskStore` (``ctx.masks``) so the
canvas — and the analyses started from this view's sidebar, which measure a set on a spectrum
into an object of their own (an abundance result, a balance measurement) — all read one
source of truth;
every mutation broadcasts ``"masks-changed"`` and this tool repaints from it.

The canvas draws the current set while the view is open: each mask's synth/visual window as a
faint red band, its fit window amber, and the active mask accented. Clicking a row snaps the
view to that mask's synth window (plus a small margin) and makes it the editable one — its
window is changed *only* by dragging the box's left/right edge handles, so an existing mask
can never be overwritten by a stray drag across the plot. A mask that does not exist yet is
drawn with the Select range tool's range gesture (:mod:`rangetool`): **+ New mask** hands over to it
in add mode, and the drawn window comes back through :meth:`MaskTool.add_window`; the same
tool starts a brand-new set from a range through :meth:`MaskTool.create_from_window`.

An edited set is flagged as modified against the file it came from; the sidebar offers a revert
to that file and a "save as" that writes a new one and adopts it as the new origin.
The toolbar's **Edit → Duplicate** (:mod:`editmenu`) adds an unsaved copy of the set to
experiment on.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

import astropy.units as u
import numpy as np
from bokeh.models import TextInput

from ..ui.controls import Field, Form
from ..ui.htmlbridge import parse_action, pulse
from tumult.linemask import is_element
from .masksources import AutoSource, BlankSource, FileSource
from ..ui.tableui import check_header
from ..ui.tool import Tool
from ..ui.units import to_range

# fraction of the synth window added on each side when a mask row is clicked
JUMP_MARGIN = 0.05

# The Auto-generate source (:class:`masksources.AutoSource`) is switched off for now: what it
# produces — one window per candidate line — has not been reconciled with the
# imported/hand-drawn view, which would want overlapping windows merged into one. Everything
# it needs is still built and wired; flipping this back on restores the pane (and the
# template's, see the `linemask_auto` flag it is passed to).
AUTO_MASKS = False


class MaskTool(Tool):
    name = "linemask"

    def __init__(self, ctx, store, input_dir, *, list_observations=None,
                 get_observation=None, master_path="", profiles=(), place_window=None,
                 on_fit_balance=None, on_measure_abund=None, on_synth=None, on_fit=None):
        super().__init__(ctx)
        self.store = store                      # shared MaskStore (also ctx.masks)
        # what the sidebar's analysis actions hand the current set to — each opens its own
        # modal on it (ionbalance.IonBalanceTool, abundancetool.AbundanceTool, the synthesis
        # modal). Injected, so the mask view stays unaware of the analyses.
        self.on_fit_balance = on_fit_balance
        self.on_measure_abund = on_measure_abund
        self.on_synth = on_synth
        self.on_fit = on_fit
        # ``place_window(set_name, callback)`` arms the Select range tool's range gesture in add mode;
        # the drawn window arrives as ``callback(xmin, xmax)`` in Å (see rangetool.RangeTool)
        self.place_window = place_window
        self.input_dir = str(input_dir)
        self.own_panel("linemask-panel-content")   # sidebar: pad / save controls
        self.include_filter = False             # the table shows only switched-on masks
        self._table_bridge = ctx.bridges.html("linemask-table")   # centre: the mask table
        # the left-sidebar list of imported mask sets (the line-mask objects)
        self._objects_bridge = ctx.bridges.html("linemask-objects")
        # the sidebar panel AND the centre bottom-panel (mask table below the spectrum, like
        # the abundance tool's #abund-bottom) toggle together. There is no toolbar button:
        # the view is entered by clicking a line-mask object in the left sidebar.
        self.own_mode_bridges((".sidebar-right", "linemask-active"),
                              ("#main", "linemask-mode"))

        # per-set settings: which element the windows diagnose, and the synth-window padding
        # (Å) that re-derives every mask's [smin,smax] from its fit window
        self.set_form = Form("linemask", [
            Field("element", "Element", kind="text", value=store.mset.element,
                  placeholder="e.g. Fe"),
            Field("pad", "Synth padding [Å]", kind="number", value=store.pad, min=0.0, step=0.1),
        ], on_change=self._on_field)

        # the ways to make a set, one pane each in the "New line mask" modal
        self.sources = [FileSource(store, self._on_source, input_dir=input_dir),
                        BlankSource(store, self._on_source)]
        if AUTO_MASKS and list_observations is not None:
            self.sources.append(AutoSource(
                store, self._on_source, ctx=ctx, list_observations=list_observations,
                get_observation=get_observation, master_path=master_path, profiles=profiles))
        self.modal = ctx.bridges.dialog("linemask-modal")   # closed once a set exists

        # per-row / save / revert actions come back on a dedicated space-delimited channel
        self.actions = TextInput(value="", visible=False, name="linemask-actions")
        self.actions.on_change("value", self._on_action)

        ctx.subscribe("masks-changed", lambda **_: self._on_masks_changed())
        # another list took (or gave back) the single left-sidebar selection: restyle our rows
        ctx.subscribe("focus-changed", lambda **_: self._render_objects())

    # --- doc wiring -----------------------------------------------------------------------
    @property
    def models(self):
        return [self.actions, *self.set_form.models,
                *(m for src in self.sources for m in src.models)]

    @property
    def roots(self):
        """Document roots the sources own (their embedded pickers)."""
        return [root for src in self.sources for root in src.roots]

    def source(self, key):
        return next((s for s in self.sources if s.key == key), None)

    # --- Tool protocol --------------------------------------------------------------------
    def activate(self):
        self.active = True
        self._push_focus("linemask")   # before the render, so the sidebar row highlights
        self.canvas.on_window_resize = self._on_edge_drag
        self._attach_annotations()
        self._render_canvas()
        self._render_panels()
        self.push_mode("on")

    def deactivate(self):
        if not self.active:
            return
        self.active = False
        if self.ctx.annot is not None:
            self.ctx.annot.detach(self.name)
        self.canvas.on_window_resize = None
        self.canvas.clear_windows()
        self.push_mode("off")
        self._push_focus("spectrum")
        self._render_objects()       # drop the sidebar highlight

    def update(self):
        """Repaint the windows after a unit or selection change."""
        auto = self.source("auto")
        if auto is not None:
            auto.refresh()       # the Auto-generate pane lists the loaded spectra
        if self.active:
            self._render_canvas()

    # --- left-sidebar line-mask objects ----------------------------------------------------
    def _push_focus(self, kind):
        """Hand the single left-sidebar selection to (or back from) the mask list."""
        if self.ctx.workspace.focus_kind == kind:
            return
        self.ctx.workspace.focus(kind)

    def select_set(self, uid):
        """A line-mask object was clicked: make it the current set and open the view."""
        try:
            uid = int(uid)
        except (TypeError, ValueError):
            return
        if self.store.find(uid) is None:
            return
        self.store.set_current(uid)
        self.ctx.tools.activate("linemask")   # stands down whatever else owned the sidebar
        self._jump_to_active()

    def remove_set(self, uid):
        try:
            uid = int(uid)
        except (TypeError, ValueError):
            return
        was_current = uid == self.store.current_uid
        placed = self.store.remove_set(uid)
        if was_current and self.active:
            self.ctx.tools.deactivate_active()
        return placed

    def objects_html(self):
        """The left-sidebar list of imported mask sets (empty string when there are none)."""
        sets = self.store.listed
        if not sets:
            return ""
        focused = self.ctx.workspace.focus_kind == "linemask"
        rows = []
        for s in sets:
            sel = " selected" if (focused and s.uid == self.store.current_uid) else ""
            name = html.escape(s.name)
            n = len(s.masks)
            # a dot marks a set edited away from the file it was imported from
            dot = '<span class="lm-dot" title="Edited since import">●</span>' if s.dirty else ""
            # the badge carries the element the set diagnoses, so the list reads at a glance
            badge = (f'<span class="result-badge maskset-badge">{html.escape(s.element)}</span>'
                     if s.element else
                     '<span class="result-badge maskset-badge maskset-badge-none" '
                     'title="No element set"><i class="ti ti-columns"></i></span>')
            rows.append(
                f'<div class="object-item maskset-item{sel}" data-maskset="{s.uid}" '
                f'title="{name} — {n} mask{"" if n == 1 else "s"}'
                f'{(" · " + html.escape(s.element)) if s.element else ""}'
                f'{" (modified)" if s.dirty else ""}">'
                f'{badge}'
                f'<span class="object-name">{name}</span>{dot}'
                f'<span class="maskset-count">{n}</span>'
                f'<button class="result-btn" data-maskset-remove="{s.uid}" '
                f'title="Remove from workspace"><i class="ti ti-x"></i></button></div>')
        return ('<div class="section-header result-section-head">Line masks</div>'
                + "".join(rows))

    # --- reactions ------------------------------------------------------------------------
    def _on_masks_changed(self):
        self.set_form.fields["pad"].value = self.store.pad
        self.set_form.fields["element"].value = self.store.mset.element
        self._render_panels()
        if self.active:
            self._attach_annotations()   # the set's element may have just changed
            self._render_canvas()

    def _attach_annotations(self):
        """Point the shared line-position overlay at this set's element, so the toolbar's
        Show/Hide button marks exactly the transitions these windows are meant to diagnose.

        The lines are confined to the set's own windows: outside them this view has nothing to
        say, and a name over a stretch of spectrum no mask covers would only read as noise. A
        set with no element yet has nothing to show at all, so it gets no button either; typing
        one into the sidebar re-runs this (through ``masks-changed``) and the button appears."""
        if self.ctx.annot is None:
            return
        element = self.store.mset.element
        if not element:
            self.ctx.annot.detach(self.name)
            return
        # the synth window is the mask's full extent on the canvas (the fit window is what it
        # falls back to when a mask carries no synth bounds, exactly as jumping does)
        windows = [(m.smin, m.smax) if np.isfinite(m.smin) and np.isfinite(m.smax)
                   else (m.lmin, m.lmax) for m in self.store.masks]
        self.ctx.annot.attach(self.name, element=element, windows=windows)

    def _render_canvas(self):
        self.canvas.highlight_windows(self.store.masks, self.store.active,
                                      editable=self.active)

    def _render_objects(self):
        self._objects_bridge.value = self.objects_html()

    def _render_panels(self):
        """Repaint the sidebar controls, the sidebar object list and the centre tables."""
        self.push_panel()
        self._render_objects()
        self._table_bridge.value = self.table_html()

    def _jump_to_active(self):
        """Snap the spectrum to the active mask's synth window (plus a 5% margin each side) so
        clicking a table row reveals the whole region the line is synthesised over."""
        m = self.store.active_mask()
        if m is None:
            return
        lo, hi = (m.smin, m.smax)
        if not (np.isfinite(lo) and np.isfinite(hi)):
            lo, hi = m.lmin, m.lmax          # fall back to the fit window (e.g. no synth bounds)
        self.canvas.focus_range(lo, hi, margin=JUMP_MARGIN)

    def _on_field(self, key, value):
        """An edit in the sidebar's per-set form (element / synth padding)."""
        if key == "element":
            value = str(value or "").strip()
            if value and not is_element(value.capitalize()):
                self.notify(f"\u201c{value}\u201d is not an element symbol.")
                self._render_panels()          # put the previous value back in the field
                return
            self.store.set_element(value)
            return
        try:
            self.store.set_pad(float(value))
        except (TypeError, ValueError):
            pass

    def _start_placing(self):
        """"+ New mask": hand over to the Select range tool's range gesture in add mode. The window
        comes back through :meth:`add_window`, which reopens this view."""
        if self.place_window is None:
            return
        self.store.set_active(None)    # nothing is editable while the new window is drawn
        self.place_window(self.store.mset.name, self.add_window)

    def add_window(self, lo, hi):
        """Append the drawn window (Å) to the open set and return to this view."""
        self.store.add(lo, hi, lam0=0.5 * (lo + hi))
        self.ctx.tools.activate("linemask")
        self.notify(f"Added mask {lo:.3f}–{hi:.3f} Å.", level="success")

    def create_from_window(self, lo, hi):
        """The Select range tool's "New line mask": a fresh set holding the drawn window (Å),
        opened in this view so its element can be set."""
        mset = self.store.add_set(BlankSource.DEFAULT_NAME)
        self.store.add(lo, hi, lam0=0.5 * (lo + hi))
        self.ctx.tools.activate("linemask")
        self._jump_to_active()
        self.notify(f"Started line mask \u201c{mset.name}\u201d with {lo:.3f}–{hi:.3f} Å.",
                    level="success")

    def _on_edge_drag(self, x0, x1):
        """An edge handle of the active mask's box was dragged (display units): resize that
        mask's fit window; its synth window follows the set padding."""
        if not self.active or self.store.active is None:
            return
        lo, hi = to_range(x0, x1, self.canvas.xunit, u.AA)
        self.store.set_fit(self.store.active, lo, hi)

    def _on_source(self, source, mset, message=""):
        """A :class:`masksources.MaskSource` produced a set (or failed): dismiss the modal,
        open the new set's view, and report. One landing for every way of adding a set, so
        each source only has to say *what* it made."""
        if mset is None:
            self.notify(message or "Could not create the line mask.", level="danger")
            return
        pulse(self.modal, "close")
        self.ctx.tools.activate("linemask")
        if mset.masks:
            self._jump_to_active()
        elif self.canvas.spectra:
            self._start_placing()
        else:
            message += " Open a spectrum to draw windows."
        self.notify(message or f"Added line mask \u201c{mset.name}\u201d.", level="success")

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "select":
            self.store.set_active(_int(arg))
            self._jump_to_active()
        elif kind == "del":
            self.store.remove(_int(arg))
        elif kind == "new":
            self._start_placing()
        elif kind == "incl":
            i = _int(arg)
            m = self.store.masks[i] if 0 <= i < len(self.store.masks) else None
            if m is not None:
                self.store.set_include([i], not m.include)
        elif kind == "check-filter":
            self.include_filter = not self.include_filter
            self._render_panels()
        elif kind in ("check-all", "check-none"):
            self.store.set_include(self._shown(), kind == "check-all")
        elif kind == "new-set":
            self._submit_source("blank", arg)
        elif kind == "auto-run":
            self._submit_source("auto", arg)
        elif kind == "auto-ref":
            auto = self.source("auto")
            if auto is not None:
                auto.set_ref(arg)
        elif kind == "fit-balance":
            self._hand_over(self.on_fit_balance)
        elif kind == "measure-abund":
            self._hand_over(self.on_measure_abund)
        elif kind == "synthesize":
            self._hand_over(self.on_synth)
        elif kind == "fit":
            self._hand_over(self.on_fit)
        elif kind == "save-as":
            self._save_as(arg)
        elif kind == "revert":
            self._revert()
        elif kind == "set-select":
            self.select_set(arg)
        elif kind == "set-remove":
            self.ctx.removal.request("linemask", arg)

    def _submit_source(self, key, arg=""):
        """A source pane's action button fired (the picker panes act on their own)."""
        source = self.source(key)
        if source is not None:
            source.submit(arg)

    def _hand_over(self, analysis):
        """Hand the current set to an analysis (which opens its own modal on it)."""
        if analysis is not None:
            analysis(self.store.mset)

    def _save_as(self, name):
        """Write the current set to ``name`` (a bare file name lands next to the source file)
        and adopt it as the set's origin, so it reads as unmodified again."""
        if not self.store.masks:
            self.notify("No line masks to save — import a file first.")
            return
        default = self.store.default_save_path()
        target = Path(str(name).strip()).expanduser() if str(name).strip() else default
        if not target.is_absolute():
            target = default.parent / target
        try:
            out = self.store.save_as(target)
        except Exception as exc:
            logging.exception("Line-mask save failed")
            self.notify(f"Could not save line masks: {exc}", level="danger")
            return
        self.notify(f"Saved {len(self.store.masks)} line masks to {out.name}.", level="success")

    def _revert(self):
        """Drop every edit and go back to the file the set was imported from / saved to."""
        if not self.store.revert():
            self.notify("This set was never imported from a file — nothing to revert to.")
            return
        self.notify(f"Reverted to {html.unescape(self.store.mset.name)}.", level="success")

    # --- rendering ------------------------------------------------------------------------
    def panel_html(self):
        """Sidebar controls: which set is open (and whether it has been edited away from its
        file), the synth padding, and the save-as / revert / clear actions. The mask table
        itself lives in the centre bottom panel, see :meth:`table_html`."""
        mset = self.store.mset
        n = len(mset.masks)
        src = mset.source_path
        dirty = mset.dirty
        badge = ('<span class="lm-dirty" title="Edited since it was imported">'
                 'Modified</span>' if dirty else "")
        srcline = (f'<div class="hint">{n} mask{"" if n == 1 else "s"} from '
                   f'{html.escape(Path(src).name)}.</div>' if src and n else
                   f'<div class="hint">{n} mask{"" if n == 1 else "s"}, not saved to a file '
                   f'yet.</div>')
        # the element is what makes a set usable downstream, so say when it is missing
        if not mset.element:
            srcline += ('<div class="hint hint-warn">Set the element these windows diagnose '
                        '— the abundance workflows need it.</div>')
        actions = ""
        if n:
            default = html.escape(self.store.default_save_path().name)
            revert = ('<button class="text-btn lm-revert">Revert</button>' if dirty else
                      '<button class="text-btn" disabled title="No edits to revert">'
                      'Revert</button>')
            actions = (
                '<div class="section-header">Save a copy</div>'
                f'<input class="field-input lm-name" type="text" value="{default}" '
                'placeholder="file name\u2026" title="Written next to the source file unless '
                'you type a full path">'
                '<div class="lm-actions">'
                '<button class="text-btn lm-saveas">Save as\u2026</button>'
                f'{revert}</div>')
        return (f'<div class="lm-title"><span class="section-header">'
                f'{html.escape(mset.name)}</span>{badge}</div>'
                + srcline + self.set_form.html() + self._analysis_html(mset) + actions)

    def _analysis_html(self, mset):
        """The actions that take this set somewhere else — one per analysis it can feed.

        A set is *the* input to the abundance workflows, so starting one belongs here rather
        than in a tool the user has to find first and then point back at this set. What comes
        back is an object of its own (a measurement, see :mod:`measurement`, or an abundance
        result), because everything it holds depends on the spectrum it was measured on, not
        on these windows.
        """
        # what each analysis does is stated in the header of the modal it opens; the flag
        # marks the analyses that diagnose one element and so need the set to name it
        actions = [(self.on_synth, "lm-synthesize", "Synthesise windows\u2026",
                    "Synthesise a spectrum over these windows", False),
                   (self.on_fit, "lm-fit", "Fit windows\u2026",
                    "Fit a target spectrum inside these windows", False),
                   (self.on_measure_abund, "lm-measure-abund", "Measure abundances",
                    "Measure these windows on a reference and invert each EW to an abundance",
                    True),
                   (self.on_fit_balance, "lm-fit-balance", "Fit ionisation balance",
                    "Measure these windows on a spectrum and fit stellar parameters", True)]
        actions = [a for a in actions if a[0] is not None]
        if not actions:
            return ""
        no_windows = ("Switch on at least one window."
                      if not any(m.include for m in mset.masks) else "")
        no_element = "" if mset.element else "Set the element these windows diagnose first."
        out = ['<div class="section-header">Analysis</div>']
        hints = []
        for _, cls, label, title, needs_element in actions:
            why = no_windows or (no_element if needs_element else "")
            out.append(f'<button class="text-btn {cls}"{" disabled" if why else ""} '
                       f'title="{html.escape(why or title)}">{label}</button>')
            if why and why not in hints:
                hints.append(why)
        out += [f'<div class="hint">{html.escape(why)}</div>' for why in hints]
        return "".join(out)

    def _shown(self):
        """Store indices of the table's rows, in wavelength order."""
        masks = self.store.masks
        order = sorted(range(len(masks)), key=lambda i: (
            masks[i].lam0 if np.isfinite(masks[i].lam0) else 0.5 * (masks[i].lmin + masks[i].lmax)))
        return [i for i in order if masks[i].include or not self.include_filter]

    def table_html(self):
        """The mask-window table: one row per mask (click to make it the editable one), a
        per-mask on/off toggle, followed by a button that starts a new mask. The button sits
        below the table's own scroll box: a sticky table footer stops a subpixel short of the
        scroll edge (its cell cannot leave the fractional-height table), so rows showed
        beneath it.

        This is exactly what a save writes out: the windows themselves. What a spectrum says
        about them belongs to the analysis that measured it, not to this set.
        """
        masks = self.store.masks
        add_row = ('<div class="lm-add-row"><button class="text-btn lm-new">'
                   '＋ New mask — then drag a range on the spectrum</button></div>')
        if not masks:
            return ('<div class="lm-scroll"><table class="linelist-table lm-table lm-table-empty">'
                    '<tbody><tr><td class="lm-empty">This set is empty — add a mask '
                    'below, then drag its edges on the spectrum to size the window.</td>'
                    f'</tr></tbody></table></div>{add_row}')
        rows = []
        # each row keeps its store index as data-idx, which is what the row actions address
        shown = self._shown()
        for i in shown:
            m = masks[i]
            sel = " lm-active" if i == self.store.active else ""
            off = "" if m.include else " lm-off"
            lam0 = f"{m.lam0:.3f}" if m.lam0 == m.lam0 else "—"   # NaN check
            inc = " checked" if m.include else ""
            rows.append(
                f'<tr class="lm-row{sel}{off}" data-idx="{i}">'
                f'<td><input type="checkbox" class="lm-incl" data-idx="{i}"{inc} '
                f'title="Use this mask (and write it when saving)"></td>'
                f'<td>{lam0}</td>'
                f'<td>{m.lmin:.3f} – {m.lmax:.3f}</td>'
                f'<td>{m.smin:.3f} – {m.smax:.3f}</td>'
                f'<td class="lm-note">{html.escape(m.note)}</td>'
                f'<td><button class="icon-btn lm-del" data-idx="{i}" title="Remove">'
                f'<i class="ti ti-x"></i></button></td></tr>')
        if not rows:
            rows.append('<tr><td class="lm-empty" colspan="6">Every mask is switched off.</td></tr>')
        head = check_header("", "include", [masks[i].include for i in shown],
                            self.include_filter, noun="masks")
        return (
            '<div class="lm-scroll"><table class="linelist-table lm-table">'
            f'<thead><tr>{head}<th>λ₀ [Å]</th>'
            '<th>fit window [Å]</th><th>synth window [Å]</th>'
            '<th>note</th><th></th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody>'
            f'</table></div>{add_row}')


def _int(arg, default=-1):
    try:
        return int(arg)
    except (TypeError, ValueError):
        return default
