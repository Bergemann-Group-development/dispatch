"""Interactive spline continuum normalisation.

Clicking the toolbar **Normalise** button enters an edit mode: the canvas shows only
the selected spectrum in its raw (un-normalised) units, with a draggable continuum
spline drawn on top. The inspector swaps to a Normalise panel that can *initialise*
the spline (auto-place knots at the flux peaks) and tune the peak spacing. Knots are
added / dragged / deleted directly on the canvas (a ``PointDrawTool``); the spectrum
is divided by the monotonic spline through them, so leaving the mode shows the
continuum-normalised spectrum (continuum at ~1).
"""
from __future__ import annotations

import logging

import numpy as np
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, PointDrawTool
from scipy.ndimage import uniform_filter1d
from scipy.signal import argrelextrema

from tumult import Spectrum
from tumult.utils import pchip_mono

from canvas import ACCENT
from controls import Field, Form
from debounce import Debouncer
from units import to_value

RED = "#e05c5c"
N_SPLINE = 400   # samples of the drawn continuum curve


class NormaliseTool:
    def __init__(self, canvas, on_panel=None, on_mode=None, notify=None, debounce_ms=200):
        self.canvas = canvas
        self.main = canvas.main
        self.on_panel = on_panel    # push the Normalise panel HTML
        self.on_mode = on_mode      # push "on"/"off" to show/hide the panel + button state
        self.notify = notify or (lambda *a, **k: None)
        self.active = False
        self._uid = None            # the object being edited
        self._suppress = False      # guard programmatic knot writes against the refit
        self._saved_range = None    # main x-range to restore when the tool closes
        self._fit_spectrum = None   # rendered spectrum the spline is fit to (smoothed+selected)
        self._debounce = Debouncer(debounce_ms)

        # `prenorm` assumes the spectrum is already normalised (continuum ≡ 1) and disables
        # the rest. `res`/`win`/`pct` only shape the rendered spectrum the spline is fit to
        # (a low-pass smooth, then a running-mean outlier rejection); `peak` (and the
        # Initialise button) is what actually places the spline knots.
        self.form = Form("norm", [
            Field("prenorm", "Already normalised", kind="checkbox", value=False),
            Field("res", "Resolution cut R=λ/Δλ (0=off)", kind="number", value=10000.0, min=0),
            Field("win", "Outlier window [Å]", kind="number", value=50.0, min=0),
            Field("pct", "Keep percentile (0=off)", kind="number", value=98.0, min=0, max=100),
            Field("peak", "Peak spacing [Å]", kind="number", value=5.0, min=0.1),
        ], on_change=self._on_field)

        # the drawn continuum curve + the draggable knots, both in display x / native flux
        self.spline_source = ColumnDataSource(dict(x=[], y=[]))
        self.spline_line = self.main.line("x", "y", source=self.spline_source,
                                          color=ACCENT, line_width=1.5, line_dash="dashed",
                                          level="overlay", visible=False)
        self.knots = ColumnDataSource(dict(x=[], y=[]))
        self.knot_glyph = self.main.scatter("x", "y", source=self.knots, size=9,
                                            fill_color=RED, line_color="#fff",
                                            line_width=1.5, level="overlay", visible=False)
        self.draw = PointDrawTool(renderers=[self.knot_glyph], add=True)
        self.main.add_tools(self.draw)
        self.knots.on_change("data", self._on_knots)

    @property
    def models(self):
        return [*self.form.models]

    # --- activation -------------------------------------------------------
    def activate(self):
        obj = self.canvas.selected()
        if obj is None:
            self.notify("Load and select a spectrum before normalising.")
            return
        self.active = True
        self._uid = obj.uid
        self._seed_knots(obj)
        self.form.fields["prenorm"].value = obj.already_normalised
        self.canvas.set_norm_edit(obj.uid)   # show this spectrum raw, hide the others
        self._zoom_full()                    # fit the whole spectrum for continuum work
        self._sync_prenorm(obj)              # render + show/hide the spline per prenorm
        if self.on_panel:
            self.on_panel(self)
        if self.on_mode:
            self.on_mode("on")

    def deactivate(self):
        if not self.active:
            return
        self.active = False
        self._uid = None
        self.spline_line.visible = False
        self.knot_glyph.visible = False
        if self.main.toolbar.active_drag is self.draw:
            self.main.toolbar.active_drag = "auto"
        self.canvas.set_norm_edit(None)      # back to the normalised view
        self._restore_range()
        if self.on_mode:
            self.on_mode("off")

    def _zoom_full(self):
        r = self.main.x_range
        self._saved_range = (r.start, r.end)
        lo, hi = self.canvas.display_bounds()
        r.start, r.end = lo, hi

    def _restore_range(self):
        if self._saved_range:
            s, e = self._saved_range
            if s is not None and e is not None:
                self.main.x_range.start, self.main.x_range.end = s, e
        self._saved_range = None

    def clear(self):
        """Drop the spline and restore the initial constant (flat) normalisation."""
        obj = self._obs()
        if obj is None:
            return
        obj.norm_knots = None
        self._set_knots([], [])
        self.canvas.normalise_object(obj)   # back to the flat percentile continuum
        self._render_selection()
        self._redraw_spline(obj)                 # show the flat initial-continuum guide
        if self.on_panel:
            self.on_panel(self)

    def _sync_prenorm(self, obj):
        """Render the view and show/hide the spline affordances for obj.already_normalised."""
        on = obj.already_normalised
        self.spline_line.visible = not on
        self.knot_glyph.visible = not on
        self.main.toolbar.active_drag = "auto" if on else self.draw
        if on:
            self.canvas.set_norm_view(None)      # show the (already-normalised) raw spectrum
        else:
            self._render_selection()
            self._redraw_spline(obj)

    def _set_prenorm(self, on):
        """Toggle 'already normalised' for the edited spectrum: continuum ≡ 1, the rest off."""
        obj = self._obs()
        if obj is None:
            return
        self.form.fields["prenorm"].value = on   # keep the panel's enable/disable state in sync
        obj.already_normalised = on
        if on:
            obj.norm_knots = None
            self._set_knots([], [])
        self.canvas.normalise_object(obj)
        self._sync_prenorm(obj)
        if self.on_panel:
            self.on_panel(self)

    def _build_fit_spectrum(self, obj):
        """The spectrum the spline is fit to and rendered: the raw spectrum low-pass
        smoothed by the resolution cut, then with outliers dropped — points whose
        excess over a running mean (taken over the `win` Å window) is above the
        `pct`-th percentile of those excesses are removed."""
        res = self.form.value("res")
        base = obj.raw.lowpass(res) if res and res > 0 else obj.raw
        win = self.form.value("win")
        pct = self.form.value("pct")
        if win and win > 0 and pct and 0 < pct < 100:
            x, y = base.values
            dx = abs(float(np.mean(np.diff(x)))) if len(x) > 1 else 0.0
            nstep = max(1, int(round(win / dx))) if dx else 1
            excess = y - uniform_filter1d(y, nstep, mode="nearest")  # vs the local mean
            mask = excess <= np.percentile(excess, pct)  # keep points within the percentile
            if mask.sum() >= 2:
                base = Spectrum(base.xdata[mask], y[mask],
                                yunit=base.yunit, xref=base.xref)
        return base

    def _render_selection(self):
        """Update the rendered spectrum (resolution cut + percentile selection) only —
        the spline knots are left untouched."""
        obj = self._obs()
        if obj is None:
            return
        self._fit_spectrum = self._build_fit_spectrum(obj)
        self.canvas.set_norm_view(self._fit_spectrum)

    def _on_field(self, key, value):
        """`prenorm` toggles the constant normalisation; `peak` re-runs the spline fit;
        `res`/`win`/`pct` only re-render the selection."""
        if not self.active or self._obs() is None:
            return
        if key == "prenorm":
            self._set_prenorm(bool(value))
        elif key == "peak":
            self._debounce("fit", self.initialise)
        else:
            self._debounce("render", self._render_selection)

    def initialise(self):
        """Place continuum knots at the flux peaks (spaced ≥ `peak` apart) of the
        current rendered selection — i.e. the resolution-cut, percentile-selected
        spectrum."""
        obj = self._obs()
        if obj is None:
            return
        self._render_selection()              # ensure the fit spectrum is current
        x, y = self._fit_spectrum.values
        dx = float(np.mean(np.diff(x))) if len(x) > 1 else 1.0
        nstep = max(1, int(abs(self.form.value("peak")) / abs(dx))) if dx else 1
        ex = argrelextrema(y, np.greater_equal, order=nstep)[0]
        ex = np.unique(np.concatenate(([0], ex, [len(x) - 1])))  # always anchor the endpoints
        kx = to_value(x[ex], self._fit_spectrum.xunit, self.canvas.xunit)  # native -> display
        self._set_knots(kx, y[ex])
        self._recompute(obj)
        self._force_repaint()
        if ex.size <= 2:  # only the endpoints survived (e.g. a strong resolution cut)
            self.notify("Few continuum peaks found — widen the peak spacing or lower the "
                        "resolution cut.", level="info")

    def _force_repaint(self):
        """Nudge the x-range a sub-pixel amount (then restore it) on the next tick. The
        WebGL knot scatter does not repaint after a programmatic empty->populate of its
        source until the next viewport change, so without this the freshly placed knots
        stay invisible until the user pans/zooms. `start` is moved *inward* so the nudge
        never collides with the range bounds (both ends sit on them when zoomed to full)."""
        def nudge():
            r = self.main.x_range
            if r.start is None or r.end is None:
                return
            s = r.start
            r.start = s + (r.end - s) * 1e-4
            curdoc().add_next_tick_callback(lambda: setattr(r, "start", s))
        curdoc().add_next_tick_callback(nudge)

    def update(self):
        """From main.update_all: leave the mode if the selection changed underneath us."""
        if self.active and (self.canvas.selected() is None
                            or self.canvas.selected().uid != self._uid):
            self.deactivate()

    # --- panel ------------------------------------------------------------
    def panel_html(self):
        if not self.active:
            return ""
        f = self.form.fields
        prenorm = bool(f["prenorm"].value)
        for k in ("res", "win", "pct", "peak"):   # everything else is moot when pre-normalised
            f[k].disabled = prenorm
        dis = " disabled" if prenorm else ""
        fields_html = (
            f["prenorm"].render("norm")
            + f["res"].render("norm")
            + '<div class="field-row">'
            + f["win"].render("norm") + f["pct"].render("norm")
            + '</div>'
            + f["peak"].render("norm")
        )
        return (
            f'{fields_html}'
            f'<button class="primary-btn norm-init"{dis}>Initialise spline</button>'
            f'<button class="text-btn norm-clear"{dis}>Clear spline</button>'
            f'<div class="hint">Resolution cut and the running-mean outlier window / '
            f'percentile shape the rendered spectrum the spline is fit to. Initialise (or '
            f'the peak spacing) places the knots; drag them to adjust, click to add, '
            f'select + Backspace to delete.</div>'
        )

    # --- knots -> spline -> normalisation ---------------------------------
    def _on_knots(self, attr, old, new):
        if self._suppress or not self.active:
            return
        self._debounce("refit", self._refit)

    def _refit(self):
        obj = self._obs()
        if obj is not None:
            self._recompute(obj)
            if self.on_panel:
                self.on_panel(self)

    def _recompute(self, obj):
        """Rebuild obj.norm_knots (native) from the on-canvas knots and renormalise."""
        d = self.knots.data
        kx_disp = np.asarray(d.get("x", []), float)
        ky = np.asarray(d.get("y", []), float)
        kx_nat = to_value(kx_disp, self.canvas.xunit, obj.raw.xunit)
        order = np.argsort(kx_nat)
        kx_nat, ky = kx_nat[order], ky[order]
        uniq = np.concatenate(([True], np.diff(kx_nat) > 0))  # pchip needs increasing x
        kx_nat, ky = kx_nat[uniq], ky[uniq]
        if kx_nat.size < 2:
            obj.norm_knots = None
        else:
            obj.norm_knots = (kx_nat, ky)
        self.canvas.normalise_object(obj)
        self._redraw_spline(obj)

    def _redraw_spline(self, obj):
        x, _ = obj.raw.values
        if obj.norm_knots is None:
            # no spline yet: show the flat initial continuum (the level the spectrum
            # is currently normalised against) as the same dashed guide
            cont = obj.continuum
            if cont is None or len(x) < 2:
                self.spline_source.data = dict(x=[], y=[])
                return
            xg = np.array([x[0], x[-1]])
            yg = np.array([cont[0], cont[-1]])
        else:
            kx_nat, ky = obj.norm_knots
            xg = np.linspace(kx_nat.min(), kx_nat.max(), N_SPLINE)
            yg = pchip_mono(kx_nat, ky)(xg)
        self.spline_source.data = dict(
            x=to_value(xg, obj.raw.xunit, self.canvas.xunit), y=yg)

    # --- helpers ----------------------------------------------------------
    def _obs(self):
        return self.canvas._find(self._uid) if self._uid is not None else None

    def _seed_knots(self, obj):
        if obj.norm_knots is None:
            self._set_knots([], [])
            return
        kx_nat, ky = obj.norm_knots
        self._set_knots(to_value(np.asarray(kx_nat, float), obj.raw.xunit,
                                 self.canvas.xunit), ky)

    def _set_knots(self, kx, ky):
        self._suppress = True
        try:
            self.knots.data = dict(x=list(map(float, kx)), y=list(map(float, ky)))
        finally:
            self._suppress = False
