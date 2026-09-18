"""Equivalent-width measurement over a drawn range.

The measurement is one of the things a range drawn with the Select range tool (:mod:`rangetool`)
can be used for, so this is a sub-panel of that tool rather than a tool of its own: the
tool hands it the range, and it measures the focused spectrum there, hatches the integrated
area on the canvas and renders its fields (method, SNR, continuum uncertainty) plus the
width and reduced width into the tool's sidebar panel. Changing a field recomputes live.
"""
from __future__ import annotations

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource

from tumult.profile import profile_names
from tumult import abundance as ab
from ..ui.controls import num_input
from ..ui.latestwork import LatestWork
from ..ui.units import to_range, to_value

ACCENT = "#4a90e2"

# Profiles the EW panel offers, in menu order — a curated subset of tumult's profile registry
# (which is the single home for the valid set). Each value is the string the tumult fitter /
# ``get_weq`` expect, so it can be passed straight through.
EW_PROFILES = tuple(p for p in ("trapz", "gauss", "voigt", "asym_gauss") if p in profile_names())

# The profiles an abundance measurement may be made with: a fitted profile (not trapz) is
# required there, because the FWHM window comes from it (see tumult.abundance).
FIT_PROFILES = ["gauss", "voigt", "asym_gauss"]

# the same profile shown under the same name in every EW selector; the fitter keeps the keys
PROFILE_LABELS = {"trapz": "Direct integration", "gauss": "Gaussian", "voigt": "Voigt",
                  "asym_gauss": "Asymmetric Gaussian"}


def profile_options(names):
    """``(value, label)`` select options for the profile ``names``."""
    return [(name, PROFILE_LABELS.get(name, name)) for name in names]


class EquivalentWidth:
    def __init__(self, ctx):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.main = ctx.canvas.main
        self.profile = "voigt"
        self.selection = None     # (x0, x1) in display units while a range is measured
        self._weq = self._rew = self._weq_err = None
        # optional error inputs: SNR (photon term) and continuum-offset % (systematic term).
        # Both off by default — the EW then shows a bare width; set either and the panel
        # shows W ± σ using the same measurement path as the abundance wizards.
        self.snr = None
        self.cont_pct = 0.0
        self._fit = None          # (profile, func, xd, yd, a, b) for redrawing
        self._work = LatestWork(lambda *args: ctx.submit(*args), "measure equivalent width",
                                lambda exc: self.notify(f"Equivalent width failed: {exc}"),
                                immediate=getattr(ctx.doc, "session_context", None) is None)

        # integrated-area hatch + its continuum/profile boundary
        self.area_source = ColumnDataSource(dict(x=[], y1=[], y2=[]))
        self.main.varea(x="x", y1="y1", y2="y2", source=self.area_source,
                        fill_color=ACCENT, fill_alpha=0.10,
                        hatch_pattern="/", hatch_color=ACCENT, hatch_alpha=0.5,
                        level="underlay")
        self.main.line("x", "y1", source=self.area_source,
                       color=ACCENT, line_width=1.2, level="underlay")

        self.main.x_range.on_change("start", lambda a, o, n: self._redraw())
        self.main.x_range.on_change("end", lambda a, o, n: self._redraw())

    # notify resolves through the context, so the footer notifier main installs later is used
    @property
    def notify(self):
        return self.ctx.notify

    # --- the spectrum measured ----------------------------------------------------------
    def target(self):
        """The focused spectrum object, or ``None`` when the left-sidebar focus is not a
        spectrum (a line mask or a result is open)."""
        if self.ctx.workspace.focus_kind != "spectrum":
            return None
        return self.canvas.spectra.selected()

    def _unavailable(self):
        """Why nothing can be measured right now ("" when it can)."""
        obj = self.target()
        if obj is None:
            return "Select a spectrum in the left sidebar to measure its equivalent width."
        if self.canvas.effective(obj).yunit != u.dimensionless_unscaled:
            return "Equivalent width needs a normalised spectrum (dimensionless flux)."
        return ""

    # --- range in / clear ---------------------------------------------------------------
    def measure(self, x0, x1):
        """Measure the focused spectrum over ``(x0, x1)`` (display units)."""
        self.selection = (min(x0, x1), max(x0, x1))
        self._compute()

    def clear(self):
        self.selection = None
        self._fit = self._weq = self._rew = self._weq_err = None
        self.area_source.data = dict(x=[], y1=[], y2=[])

    def update(self):
        """The canvas changed (focus, units, filter state): re-measure the current range."""
        if self.selection is not None:
            self._compute()

    def set_profile(self, name):
        if name not in EW_PROFILES:
            return
        self.profile = name
        self.update()

    def set_snr(self, raw):
        """SNR for the photon EW error (blank clears it — no photon term)."""
        raw = str(raw).strip()
        try:
            self.snr = float(raw) if raw else None
        except ValueError:
            return
        if self.snr is not None and self.snr <= 0:
            self.snr = None
        self.update()

    def set_cont(self, raw):
        """Continuum-offset uncertainty [%] for the systematic EW error (blank / 0 clears it)."""
        raw = str(raw).strip()
        try:
            self.cont_pct = float(raw) if raw else 0.0
        except ValueError:
            return
        self.update()

    def panel_html(self):
        why = self._unavailable()
        if why:
            return f'<div class="hint">{why}</div>'
        opts = "".join(
            f'<option value="{p}"{" selected" if p == self.profile else ""}>{label}</option>'
            for p, label in profile_options(EW_PROFILES)
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
            + num_input('class="field-input ew-snr-input"', snr, min=0, placeholder="—")
            + '</div>'
            f'<div class="field-group"><label>Continuum uncertainty [%]</label>'
            + num_input('class="field-input ew-cont-input"', cont, min=0, step=0.1)
            + '</div>'
            f'<div class="field-group"><label>Equivalent width</label>'
            f'<div class="field-static ew-value">{ew}</div></div>'
            f'<div class="field-group"><label>Reduced EW (log W/λ)</label>'
            f'<div class="field-static ew-value">{rew}</div></div>'
        )

    # --- computation + drawing --------------------------------------------
    def _compute(self):
        obj = self.target()
        if obj is None or self.selection is None or self._unavailable():
            self._weq = self._rew = self._weq_err = self._fit = None
            self.area_source.data = dict(x=[], y1=[], y2=[])
            return
        spec = self.canvas.effective(obj)  # filtered or unfiltered, per the canvas state
        x0, x1 = self.selection
        native = spec.xunit
        display_xunit, display_yunit = self.canvas.xunit, self.canvas.yunit
        profile, snr, cont_pct = self.profile, self.snr, self.cont_pct
        selection = self.selection

        def compute(token):
            b0, b1 = to_range(x0, x1, display_xunit, native)
            sub = spec.crop(xmin=b0, xmax=b1)
            if len(sub) < 3:
                return "narrow", None
            token.check()
            prof = ab.measure_ew(sub, snr=snr, cont_pct=cont_pct,
                                 profile=profile, x0=0.5 * (b0 + b1))
            disp = sub.convert_units(xunit=display_xunit, yunit=display_yunit)
            xd, yd = disp.values
            if xd[0] > xd[-1]:
                xd, yd = xd[::-1], yd[::-1]
            a, b = min(x0, x1), max(x0, x1)
            fit = dict(profile=profile, func=prof.func if prof.fitted else None,
                       native=native, display_xunit=display_xunit,
                       xd=xd, yd=yd, a=a, b=b)
            return "ok", (prof.maa, prof.width_err_maa, prof.rew, fit)

        def apply(outcome):
            if (self.selection != selection or self.profile != profile or self.snr != snr
                    or self.cont_pct != cont_pct or self.target() is not obj):
                return
            status, value = outcome
            if status == "narrow":
                self.notify("Selected range is too narrow to measure an equivalent width.")
                return
            self._weq, self._weq_err, self._rew, self._fit = value
            self._redraw()
            tool = self.ctx.tools._tools.get("range")
            if tool is not None and tool.active:
                tool.push_panel()

        self._work.request(compute, apply)

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
            if f["func"] is not None:
                native_x = to_value(xg, f["display_xunit"], f["native"])
                y1 = f["func"](native_x)
            else:
                y1 = np.interp(xg, f["xd"], f["yd"])
        self.area_source.data = dict(x=xg, y1=np.asarray(y1, dtype=float),
                                     y2=np.ones_like(xg))
