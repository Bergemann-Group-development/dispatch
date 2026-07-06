"""Equivalent-width measurement.

Workflow: drag an x-range with the select tool (marks a region, no compute) →
click the EW button in the app toolbar → the inspector switches to an EW panel
showing the method selector, the equivalent width and the reduced EW, and the
integrated area is hatched on the canvas. Changing the method recomputes live.
Esc clears the selection and restores the inspector.
"""
from __future__ import annotations

import logging

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource

from tumult.profile import profile_names
from tumult.utils import fit_line_profile

import abundance as ab
from rangeselect import RangeSelect
from tool import Tool
from units import to_range

ACCENT = "#4a90e2"

# Profiles the EW tool offers, in menu order — a curated subset of tumult's profile registry
# (which is the single home for the valid set). Each value is the string the tumult fitter /
# ``get_weq`` expect, so it can be passed straight through.
EW_PROFILES = tuple(p for p in ("trapz", "gauss", "voigt", "asym_gauss") if p in profile_names())


class EquivalentWidth(Tool, RangeSelect):
    name = "ew"

    def __init__(self, ctx):
        Tool.__init__(self, ctx)
        RangeSelect.__init__(self, ctx.canvas, color=ACCENT, description="Select x-range",
                             fill_alpha=0.06)
        self.own_panel("ew-panel-content")
        self.own_mode_bridges((".sidebar-right", "ew-active"), ("#ew-btn", "active"))
        # the EW follows whichever spectrum (filtered/unfiltered) the canvas shows
        ctx.subscribe("filter-state", lambda **_: self.update())
        self.profile = "voigt"
        self._active_uid = None
        self._weq = self._rew = self._weq_err = None
        # optional error inputs: SNR (photon term) and continuum-offset % (systematic term).
        # Both off by default — the EW then shows a bare width, exactly as before; set either and
        # the panel shows W ± σ using the same measurement path as the abundance wizards.
        self.snr = None
        self.cont_pct = 0.0
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
        else:
            self.push_panel()          # show the (empty) panel until a range is drawn
        self.push_mode("on")

    def deactivate(self):
        """Leave EW mode: drop the selection, clear the hatch, deselect the tool."""
        was = self.active or self.selection is not None
        self.active = False
        self._fit = self._weq = self._rew = self._weq_err = None
        self.clear_selection()                      # hide the band, drop the range
        self.area_source.data = dict(x=[], y1=[], y2=[])
        self.disarm()                               # restore panning
        if was:
            self.push_mode("off")

    # the Esc-key path still calls cancel(); it is the same as leaving the mode
    cancel = deactivate

    def set_profile(self, name):
        if name not in EW_PROFILES:
            return
        self.profile = name
        if self.active:
            self._compute()

    def set_snr(self, raw):
        """SNR for the photon EW error (blank clears it — no photon term)."""
        raw = str(raw).strip()
        try:
            self.snr = float(raw) if raw else None
        except ValueError:
            return
        if self.snr is not None and self.snr <= 0:
            self.snr = None
        if self.active:
            self._compute()

    def set_cont(self, raw):
        """Continuum-offset uncertainty [%] for the systematic EW error (blank / 0 clears it)."""
        raw = str(raw).strip()
        try:
            self.cont_pct = float(raw) if raw else 0.0
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
            for p in EW_PROFILES
        )
        snr = "" if self.snr is None else f"{self.snr:g}"
        cont = f"{self.cont_pct:g}"
        # always show the error when we have one: fitted profiles carry a fit-covariance error
        # even with no SNR / continuum term set, and SNR / continuum-% add to it in quadrature.
        if self._weq is None:
            ew = "—"
        elif self._weq_err is not None and np.isfinite(self._weq_err):
            ew = f"{self._weq:.2f} ± {self._weq_err:.2f} mÅ"
        else:
            ew = f"{self._weq:.2f} mÅ"
        rew = (f"{self._rew:.3f}" if self._rew is not None and np.isfinite(self._rew)
               else "—")
        return (
            f'<div class="field-group"><label>Method</label>'
            f'<select class="field-input ew-method-select">{opts}</select></div>'
            f'<div class="field-group"><label>SNR</label>'
            f'<input type="number" class="field-input ew-snr-input" min="0" '
            f'placeholder="—" value="{snr}"></div>'
            f'<div class="field-group"><label>Continuum uncertainty [%]</label>'
            f'<input type="number" class="field-input ew-cont-input" min="0" step="0.1" '
            f'value="{cont}"></div>'
            f'<div class="field-group"><label>Equivalent width</label>'
            f'<div class="field-static ew-value">{ew}</div></div>'
            f'<div class="field-group"><label>Reduced EW (log W/λ)</label>'
            f'<div class="field-static ew-value">{rew}</div></div>'
        )

    # --- computation + drawing --------------------------------------------
    def _compute(self):
        obj = self.canvas.selected()
        if obj is None or self.selection is None:
            return
        spec = self.canvas.effective(obj)  # filtered or unfiltered, per the canvas state
        if spec.yunit != u.dimensionless_unscaled:  # EW needs normalized flux
            self._weq = self._rew = self._weq_err = self._fit = None
            self.area_source.data = dict(x=[], y1=[], y2=[])
            self.push_panel()
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
            # measure the width, reduced EW and (when SNR / continuum-% are set) the error split
            # through the same primitive the abundance wizards use, so the numbers agree exactly.
            prof = ab.measure_ew(sub, snr=self.snr, cont_pct=self.cont_pct,
                                 profile=self.profile, x0=0.5 * (b0 + b1))
            self._weq = prof.mAA
            self._weq_err = prof.width_err_mAA
            self._rew = prof.rew

            disp = sub.convert_units(xunit=self.canvas.xunit, yunit=self.canvas.yunit)
            xd, yd = disp.values
            if xd[0] > xd[-1]:
                xd, yd = xd[::-1], yd[::-1]
            a, b = min(x0, x1), max(x0, x1)
            func = None
            if self.profile != "trapz":
                _, _, func = fit_line_profile(xd, yd, x0=0.5 * (a + b), profile=self.profile)
            self._fit = dict(profile=self.profile, func=func, xd=xd, yd=yd, a=a, b=b)
            self._redraw()
            self.push_panel()
        except Exception:
            logging.exception("Equivalent width computation failed")

    def _redraw(self):
        """trapz hatches only the selection; fitted profiles span the whole canvas."""
        f = self._fit
        if not f:
            return
        if f["profile"] == "trapz":
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
