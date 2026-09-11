"""The Mask tool: one range gesture on the canvas, and the things a range can start.

Clicking **Mask** in the toolbar makes width-only box-select the active drag and swaps the
right sidebar to this tool's panel. Once a range is drawn its band stays shaded and the
panel shows it as editable λ min / λ max fields (Å) with the actions that take it
somewhere:

* **New line mask** — a fresh line-mask set holding this window (never an existing set).
* **Synthesise…** / **Fit…** — the synthesis modal (:mod:`synthmodal`) on this window.
* **Equivalent width** — a toggle: measure the focused spectrum over it (:mod:`ew`) and
  keep the EW fields in the panel, following the range; press again to hide them.

The LineMasks view reuses the same gesture to extend the open set: ``begin_add`` puts the
tool in *add* mode, where the only action is "Add window to <set>" and leaving the tool
returns to that view. See ``docs/GUIWorkflow.md``.
"""
from __future__ import annotations

import html

import astropy.units as u
import numpy as np

from controls import Field, Form
from ew import ACCENT, EquivalentWidth
from rangeselect import RangeSelect
from tool import Tool
from units import to_range


class RangeTool(Tool):
    name = "range"

    def __init__(self, ctx, *, on_new_mask=None, on_synth=None, on_fit=None):
        super().__init__(ctx)
        self.on_new_mask = on_new_mask   # (xmin, xmax) Å -> a new line-mask object
        self.on_synth = on_synth         # (xmin, xmax) Å -> synthesis modal
        self.on_fit = on_fit             # (xmin, xmax) Å -> fit modal
        self.range = RangeSelect(self.canvas, ACCENT, "Select a wavelength range",
                                 fill_alpha=0.06)
        self.range.on_select = self._on_select
        self.ew = EquivalentWidth(ctx)
        self.selection = None            # (xmin, xmax) in Å once a range is drawn
        self._measuring = False          # True once "Equivalent width" was pressed
        self._add = None                 # (set name, callback(xmin, xmax)) while adding
        self.own_panel("range-panel-content")
        self.own_mode_bridges((".sidebar-right", "range-active"), ("#mask-btn", "active"))
        self.window_form = Form("rangesel", [
            Field("xmin", "λ min [Å]", kind="number", value=None, min=0, inline=True),
            Field("xmax", "λ max [Å]", kind="number", value=None, min=0, inline=True),
        ], on_change=self._on_field)
        # the EW follows whichever spectrum (filtered/unfiltered) the canvas shows
        ctx.subscribe("filter-state", lambda **_: self.update())

    @property
    def models(self):
        return list(self.window_form.models)

    # --- Tool protocol ------------------------------------------------------------------
    def activate(self):
        self.active = True
        self.range.arm()
        self.push_panel()
        self.push_mode("on")

    def deactivate(self):
        if not self.active:
            return
        self.active = False
        self._add = None
        self.selection = None
        self._measuring = False
        self.range.clear_selection()
        self.range.disarm()
        self.ew.clear()
        self.push_mode("off")

    def cancel(self):
        """The panel's × or Esc: leave the tool; in add mode, back to the LineMasks view."""
        if self._add is not None:
            self.ctx.tools.activate("linemask")
        else:
            self.deactivate()

    def update(self):
        """The canvas changed (focus, units, filter): keep the band and the EW in step."""
        if not self.active or self.selection is None:
            return
        self._show_band()
        if self._measuring:
            self.ew.measure(*self._display_range())
        self.push_panel()

    def select_range(self, lo, hi):
        """Open the tool on a given range (Å), e.g. a spectrum's full extent."""
        self.ctx.tools.activate("range")
        self._set_selection(lo, hi)

    # --- add mode (the LineMasks view's "+ New mask") -----------------------------------
    def begin_add(self, set_name, callback):
        """Arm the gesture to draw one window for the open set; ``callback(xmin, xmax)``
        (Å) receives it, and is expected to reopen the LineMasks view."""
        self.ctx.tools.activate("range")
        self._add = (set_name, callback)
        self.push_panel()

    # --- the range ----------------------------------------------------------------------
    def _display_range(self):
        lo, hi = to_range(self.selection[0], self.selection[1], u.AA, self.canvas.xunit)
        return lo, hi

    def _show_band(self):
        self.range.show_band(*self._display_range())

    def _on_select(self, x0, x1):
        """A range was drawn (display units): keep it in Å, shade it, measure it."""
        if not (np.isfinite(x0) and np.isfinite(x1)) or x0 == x1:
            return
        lo, hi = to_range(x0, x1, self.canvas.xunit, u.AA)
        self._set_selection(lo, hi)

    def _on_field(self, key, value):
        """λ min / λ max edited by hand."""
        v = self.window_form.values()
        try:
            lo, hi = float(v["xmin"]), float(v["xmax"])
        except (TypeError, ValueError):
            return
        if not hi > lo > 0:
            return
        self._set_selection(lo, hi)

    def _set_selection(self, lo, hi):
        self.selection = (round(lo, 3), round(hi, 3))
        self.window_form.fields["xmin"].value = self.selection[0]
        self.window_form.fields["xmax"].value = self.selection[1]
        self._show_band()
        if self._measuring:
            self.ew.measure(*self._display_range())
        self.push_panel()

    # --- actions ------------------------------------------------------------------------
    def new_mask(self):
        if self.selection is not None and self.on_new_mask is not None:
            self.on_new_mask(*self.selection)

    def synthesise(self):
        if self.selection is not None and self.on_synth is not None:
            self.on_synth(*self.selection)

    def fit(self):
        if self.selection is not None and self.on_fit is not None:
            self.on_fit(*self.selection)

    def toggle_ew(self):
        """"Equivalent width": show the measurement of the focused spectrum over the range
        (it then follows the range), or hide it again."""
        if self.selection is None:
            return
        self._measuring = not self._measuring
        if self._measuring:
            self.ew.measure(*self._display_range())
        else:
            self.ew.clear()
        self.push_panel()

    def ew_setting(self, setter, value):
        """An EW field changed: apply it (which re-measures) and repaint the panel."""
        setter(value)
        self.push_panel()

    def add_window(self):
        if self.selection is None or self._add is None:
            return
        _, callback = self._add
        self._add = None
        callback(*self.selection)

    # --- rendering ----------------------------------------------------------------------
    def panel_html(self):
        if not self.active:
            return ""
        if self.selection is None:
            if self._add is not None:
                return (f'<div class="hint">Drag a range on the spectrum to add a window to '
                        f'<b>{html.escape(self._add[0])}</b>.</div>')
            return '<div class="hint">Drag a range on the spectrum.</div>'
        out = ['<div class="section-header">Range</div>',
               f'<div id="{self.window_form.slot_id}">{self.window_form.html()}</div>']
        if self._add is not None:
            out.append('<button class="text-btn range-add">Add window to '
                       f'{html.escape(self._add[0])}</button>')
            return "".join(out)
        out.append('<div class="section-header">Use the range</div>'
                   '<div class="range-actions">'
                   '<button class="text-btn range-new-mask" title="Start a line-mask set '
                   'holding this window">New line mask</button>'
                   '<button class="text-btn range-synth" title="Synthesise a spectrum over '
                   'this window">Synthesise…</button>'
                   '<button class="text-btn range-fit" title="Fit a spectrum over this '
                   'window">Fit…</button>'
                   f'<button class="text-btn range-ew{" active" if self._measuring else ""}" '
                   'title="Show / hide the equivalent width of the focused spectrum over this '
                   'window">Equivalent width</button></div>')
        if self._measuring:
            out.append('<div class="section-header">Equivalent width</div>')
            out.append(self.ew.panel_html())
        return "".join(out)
