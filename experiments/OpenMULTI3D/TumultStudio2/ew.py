"""Equivalent-width measurement.

Workflow: drag an x-range with the select tool (marks a region, no compute) →
click the EW button in the app toolbar → the inspector switches to an EW panel
showing the method selector, the equivalent width and the reduced EW, and the
integrated area is hatched on the canvas. Changing the method recomputes live.
Esc clears the selection and restores the inspector.
"""
from __future__ import annotations

import logging
from enum import StrEnum

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource

from tumult.utils import fit_line_profile

from rangeselect import RangeSelect
from units import to_value, to_range

ACCENT = "#4a90e2"


class Profile(StrEnum):
    """EW line-profile model. The value is the name the tumult fitter expects, so a
    member can be passed straight to ``fit_line_profile`` / ``get_weq`` as a string."""
    TRAPZ = "trapz"
    GAUSS = "gauss"
    VOIGT = "voigt"
    ASYM_GAUSS = "asym_gauss"


class EquivalentWidth(RangeSelect):
    def __init__(self, canvas, on_panel=None, on_mode=None, notify=None):
        super().__init__(canvas, color=ACCENT, description="Select x-range", fill_alpha=0.06)
        self.profile = Profile.VOIGT
        self.on_panel = on_panel  # push the EW panel HTML (method + values)
        self.on_mode = on_mode    # push "on"/"off" to show/hide the EW panel
        self.notify = notify or (lambda *a, **k: None)  # footer info message
        self.active = False
        self._active_uid = None
        self._weq = self._ewr = None
        self._fit = None          # (profile, func, xd, yd, a, b) for redrawing

        # integrated-area hatch + its continuum/profile boundary (EW mode only)
        self.area_source = ColumnDataSource(dict(x=[], y1=[], y2=[]))
        self.main.varea(x="x", y1="y1", y2="y2", source=self.area_source,
                        fill_color=ACCENT, fill_alpha=0.10,
                        hatch_pattern="/", hatch_color=ACCENT, hatch_alpha=0.5,
                        level="underlay")
        self.main.line("x", "y1", source=self.area_source,
                       color=ACCENT, line_width=1.2, level="underlay")

        self.main.x_range.on_change("start", lambda a, o, n: self._redraw())
        self.main.x_range.on_change("end", lambda a, o, n: self._redraw())

    # --- x-range selection (no compute) -----------------------------------
    def on_select(self, x0, x1):
        self.show_band(x0, x1)
        if self.active:  # already in EW mode -> recompute for the new region
            self._compute()

    # --- EW mode: activate / deactivate / method --------------------------
    def activate(self):
        """Enter EW mode: arm the x-range select tool and show the panel."""
        self.active = True
        self.arm()                                  # select the x-range tool
        obj = self.canvas.selected()
        self._active_uid = obj.uid if obj else None
        if self.selection is not None:
            self._compute()
        elif self.on_panel:
            self.on_panel(self)        # show the (empty) panel until a range is drawn
        if self.on_mode:
            self.on_mode("on")

    def deactivate(self):
        """Leave EW mode: drop the selection, clear the hatch, deselect the tool."""
        was = self.active or self.selection is not None
        self.active = False
        self._fit = self._weq = self._ewr = None
        self.clear_selection()                      # hide the band, drop the range
        self.area_source.data = dict(x=[], y1=[], y2=[])
        self.disarm()                               # restore panning
        if was and self.on_mode:
            self.on_mode("off")

    # the Esc-key path still calls cancel(); it is the same as leaving the mode
    cancel = deactivate

    def set_profile(self, name):
        try:
            self.profile = Profile(name)
        except ValueError:
            return
        if self.active:
            self._compute()

    def update(self):
        """Re-evaluate the EW for the current canvas state (selection target, and the
        filtered/unfiltered spectrum); leave EW mode if the target is no longer valid."""
        if not self.active:
            return
        obj = self.canvas.selected()
        if (obj is None or obj.uid != self._active_uid
                or self.canvas.effective(obj).yunit != u.dimensionless_unscaled):
            self.cancel()
        elif self.selection is not None:
            self._compute()

    def panel_html(self):
        if not self.active:
            return ""
        opts = "".join(
            f'<option value="{p}"{" selected" if p == self.profile else ""}>{p}</option>'
            for p in Profile
        )
        ew = f"{self._weq:.2f} mÅ" if self._weq is not None else "—"
        ewr = (f"{self._ewr:.3f}" if self._ewr is not None and np.isfinite(self._ewr)
               else "—")
        return (
            f'<div class="field-group"><label>Method</label>'
            f'<select class="field-input ew-method-select">{opts}</select></div>'
            f'<div class="field-group"><label>Equivalent width</label>'
            f'<div class="field-static ew-value">{ew}</div></div>'
            f'<div class="field-group"><label>Reduced EW (log W/λ)</label>'
            f'<div class="field-static ew-value">{ewr}</div></div>'
        )

    # --- computation + drawing --------------------------------------------
    def _compute(self):
        obj = self.canvas.selected()
        if obj is None or self.selection is None:
            return
        spec = self.canvas.effective(obj)  # filtered or unfiltered, per the canvas state
        if spec.yunit != u.dimensionless_unscaled:  # EW needs normalized flux
            self._weq = self._ewr = self._fit = None
            self.area_source.data = dict(x=[], y1=[], y2=[])
            if self.on_panel:
                self.on_panel(self)
            self.notify("Equivalent width needs a normalized spectrum (dimensionless flux).")
            return
        try:
            x0, x1 = self.selection
            native = spec.xunit
            b0, b1 = to_range(x0, x1, self.canvas.xunit, native)
            sub = spec.crop(xmin=b0, xmax=b1)
            if len(sub) < 3:
                self.notify("Selected range is too narrow to measure an equivalent width.")
                return
            weq = sub.get_weq(profile=self.profile)
            self._weq = weq.to(u.mAA).value
            lam_aa = to_value(0.5 * (b0 + b1), native, u.AA)
            wa = weq.to(u.AA).value
            self._ewr = float(np.log10(wa / lam_aa)) if wa > 0 and lam_aa > 0 else np.nan

            disp = sub.convert_units(xunit=self.canvas.xunit, yunit=self.canvas.yunit)
            xd, yd = disp.values
            if xd[0] > xd[-1]:
                xd, yd = xd[::-1], yd[::-1]
            a, b = min(x0, x1), max(x0, x1)
            func = None
            if self.profile != Profile.TRAPZ:
                _, _, func = fit_line_profile(xd, yd, x0=0.5 * (a + b), profile=self.profile)
            self._fit = dict(profile=self.profile, func=func, xd=xd, yd=yd, a=a, b=b)
            self._redraw()
            if self.on_panel:
                self.on_panel(self)
        except Exception:
            logging.exception("Equivalent width computation failed")

    def _redraw(self):
        """trapz hatches only the selection; fitted profiles span the whole canvas."""
        f = self._fit
        if not f:
            return
        if f["profile"] == Profile.TRAPZ:
            xg = np.linspace(f["a"], f["b"], 256)
            y1 = np.interp(xg, f["xd"], f["yd"])
        else:
            r = self.main.x_range
            if r.start is None or r.end is None:
                return
            lo, hi = min(r.start, r.end), max(r.start, r.end)
            xg = np.linspace(lo, hi, 512)
            y1 = f["func"](xg) if f["func"] is not None else np.interp(xg, f["xd"], f["yd"])
        self.area_source.data = dict(x=xg, y1=np.asarray(y1, dtype=float),
                                     y2=np.ones_like(xg))
