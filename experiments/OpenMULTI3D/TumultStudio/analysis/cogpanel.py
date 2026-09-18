"""The panel beside the spectrum that shows where a line's abundance came from.

Each method leaves a different curve behind, and each gets its own panel over the same
abundance axis, so the solution is seen rather than taken on trust:

:class:`CogPanel` — the curve-of-growth inversion (:func:`tumult.abundance.adaptive_cog`)
in the plane it is read in, reduced EW against abundance: the rungs the line was bracketed
with, the PCHIP interpolant :func:`tumult.abundance.cog_invert` reads the abundance off,
the measured REW with its ±1σ band and the inverted abundance with its propagated ±σ band.

:class:`Chi2Panel` — the profile fit's own objective (:func:`tumult.abundance.chi2_profile`)
against abundance, with the solution and the chi-square the table reports marked on it. A
profile fit inverts no width, so it has no curve of growth; what it minimised is the curve
worth seeing.

The abundance tool and the NLTE-correction tool both show these, over an axis that is the
element's [X/Fe] or A(X) in one and a difference ΔA in the other, so the axis label and the
title's value label are the caller's.
"""
from __future__ import annotations

import numpy as np
from bokeh.models import BoxAnnotation, ColumnDataSource, HoverTool, Range1d, Span
from bokeh.plotting import figure

from tumult.abundance import reduced_ew
from tumult.utils import pchip_mono

from ..spectrum.canvas import ACCENT, SYNTH_FIT, style_dark

CURVE_POINTS = 80            # samples of the PCHIP interpolant drawn through the rungs
SOLUTION_COLOR = "#e05c5c"   # the inverted abundance and its ±σ band


class _Panel:
    """A styled figure over an abundance axis. ``crop`` fixes the ranges to the drawn values
    (:meth:`diagnostics.DiagnosticPanels.crop`) and ``fmt`` formats a value ± error
    (:func:`resultview.format_value`)."""

    y_label = ""

    def __init__(self, name, *, crop, fmt, x_label="ΔA [dex]"):
        self.crop, self.fmt = crop, fmt
        fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                     tools="pan,wheel_zoom,box_zoom,reset", active_scroll="wheel_zoom",
                     min_border_left=62, min_border_bottom=42,
                     x_range=Range1d(-1.0, 1.0), y_range=Range1d(0.0, 1.0))
        style_dark(fig)
        fig.title.text_color = "#cfd2d6"
        fig.title.text_font_size = "11px"
        fig.name = name
        fig.xaxis.axis_label = x_label
        fig.yaxis.axis_label = self.y_label
        self.fig = fig
        self._build(fig)

    def set_x_label(self, label):
        self.fig.xaxis.axis_label = label


class CogPanel(_Panel):
    """The curve-of-growth inversion: rungs, the interpolant, the measured REW and the
    abundance read off it."""

    y_label = "Reduced EW  log(W/λ)"

    def _build(self, fig):
        self.curve = ColumnDataSource(dict(a=[], rew=[]))
        self.rungs = ColumnDataSource(dict(a=[], rew=[], ew=[]))
        fig.line("a", "rew", source=self.curve, line_color=SYNTH_FIT, line_width=1.5)
        rungs = fig.scatter("a", "rew", source=self.rungs, size=8, fill_color=SYNTH_FIT,
                            line_color="#fff", line_width=1)
        fig.add_tools(HoverTool(renderers=[rungs],
                                tooltips=[("A", "@a{+0.00} dex"), ("REW", "@rew{0.000}"),
                                          ("EW", "@ew{0.0} mÅ")]))
        # the measurement (horizontal) and its solution (vertical), each with a ±1σ band
        self.ew_band = BoxAnnotation(fill_color=ACCENT, fill_alpha=0.12, visible=False)
        self.ew_line = Span(location=0, dimension="width", line_color=ACCENT,
                            line_width=1.2, line_dash="dashed", visible=False)
        self.a_band = BoxAnnotation(fill_color=SOLUTION_COLOR, fill_alpha=0.12, visible=False)
        self.a_line = Span(location=0, dimension="height", line_color=SOLUTION_COLOR,
                           line_width=1.2, line_dash="dashed", visible=False)
        self._annotations = (self.ew_band, self.ew_line, self.a_band, self.a_line)
        for layout in self._annotations:
            fig.add_layout(layout)

    def render(self, ln, *, value_label="ΔA"):
        """Draw ``ln``'s curve of growth, or clear the panel when the line has none — it
        was not solved by this method, or was dropped before its rungs were measured."""
        a_grid = np.asarray(ln.a_grid, dtype=float)
        ew_grid = np.asarray(ln.ew_grid, dtype=float)
        if a_grid.size < 2 or a_grid.size != ew_grid.size:
            self.clear()
            return
        if not (np.all(ew_grid > 0) and ln.ew > 0):
            self.clear()               # a non-positive width has no reduced EW
            return
        rew_grid = np.array([reduced_ew(w, ln.lam0) for w in ew_grid])
        self.rungs.data = dict(a=a_grid.tolist(), rew=rew_grid.tolist(),
                               ew=ew_grid.tolist())
        # the inversion reads A(REW) off a PCHIP interpolant through the rungs sorted by
        # REW (:func:`tumult.abundance.cog_invert`); the same interpolant is drawn, sampled
        # in REW, so the curve shows exactly what the solution was read from
        order = np.argsort(rew_grid)
        x, a = rew_grid[order], a_grid[order]
        curve_rew = np.linspace(x[0], x[-1], CURVE_POINTS)
        curve_a = pchip_mono(x, a)(curve_rew)
        self.curve.data = dict(a=np.asarray(curve_a).tolist(), rew=curve_rew.tolist())
        ys = list(rew_grid)
        # the measured REW and its error, σ_REW = σ_W / (W ln10)
        rew = reduced_ew(ln.ew, ln.lam0)
        err = float(ln.ew_err) / (ln.ew * np.log(10.0)) if np.isfinite(ln.ew_err) else 0.0
        self.ew_line.visible = self.ew_band.visible = True
        self.ew_line.location = rew
        self.ew_band.bottom, self.ew_band.top = rew - err, rew + err
        ys += [rew - err, rew + err]
        a_ok = np.isfinite(ln.a)
        self.a_line.visible = self.a_band.visible = bool(a_ok)
        if a_ok:
            sig = float(ln.sigma_a) if np.isfinite(ln.sigma_a) else 0.0
            self.a_line.location = float(ln.a)
            self.a_band.left, self.a_band.right = ln.a - sig, ln.a + sig
        title = f"{ln.ion} {ln.lam0:.3f} Å"
        if a_ok:
            title += f"   {value_label} = {self.fmt(ln.a, ln.sigma_a, plus=True)}"
        self.fig.title.text = title
        self.crop(list(a_grid), ys, self.fig)

    def clear(self):
        self.rungs.data = dict(a=[], rew=[], ew=[])
        self.curve.data = dict(a=[], rew=[])
        self.fig.title.text = ""
        for layout in self._annotations:
            layout.visible = False


class Chi2Panel(_Panel):
    """The profile fit's objective against abundance, sliced through its own solution.

    Everything but the abundance is held at the values the fit stopped on, so the minimum
    drawn here is the line's abundance and the value at it is the chi-square the table
    reports (:func:`tumult.abundance.reduced_chi2`). The scan reads cached rungs, so it
    costs no synthesis — but it is one interpolation per sample, so render on demand.
    """

    y_label = "Reduced χ²"

    def _build(self, fig):
        self.curve = ColumnDataSource(dict(a=[], chi2=[]))
        fig.line("a", "chi2", source=self.curve, line_color=SYNTH_FIT, line_width=1.5)
        self.best = ColumnDataSource(dict(a=[], chi2=[]))
        best = fig.scatter("a", "chi2", source=self.best, size=9,
                           fill_color=SOLUTION_COLOR, line_color="#fff", line_width=1)
        fig.add_tools(HoverTool(renderers=[best],
                                tooltips=[("A", "@a{+0.00} dex"), ("χ²ᵣ", "@chi2{0.000}")]))
        # the solution and its ±σ band, drawn as the CoG panel draws its own
        self.a_band = BoxAnnotation(fill_color=SOLUTION_COLOR, fill_alpha=0.12, visible=False)
        self.a_line = Span(location=0, dimension="height", line_color=SOLUTION_COLOR,
                           line_width=1.2, line_dash="dashed", visible=False)
        # the chi-square reached there: the number the table carries for this line
        self.chi2_line = Span(location=0, dimension="width", line_color=ACCENT,
                              line_width=1.2, line_dash="dashed", visible=False)
        self._annotations = (self.a_band, self.a_line, self.chi2_line)
        for layout in self._annotations:
            fig.add_layout(layout)

    def render(self, ln, a_grid, chi2, *, value_label="A"):
        """Draw the scan ``(a_grid, chi2)`` of :func:`tumult.abundance.chi2_profile`, or
        clear the panel when the line carries no fit to scan."""
        a_grid = np.asarray(a_grid, dtype=float)
        chi2 = np.asarray(chi2, dtype=float)
        finite = np.isfinite(a_grid) & np.isfinite(chi2)
        if finite.sum() < 2:
            self.clear()
            return
        a_grid, chi2 = a_grid[finite], chi2[finite]
        self.curve.data = dict(a=a_grid.tolist(), chi2=chi2.tolist())
        a_ok = np.isfinite(ln.a)
        self.a_line.visible = self.a_band.visible = bool(a_ok)
        best = float("nan")
        if a_ok:
            sig = float(ln.sigma_a) if np.isfinite(ln.sigma_a) else 0.0
            self.a_line.location = float(ln.a)
            self.a_band.left, self.a_band.right = ln.a - sig, ln.a + sig
            best = float(np.interp(ln.a, a_grid, chi2))
        self.best.data = (dict(a=[float(ln.a)], chi2=[best]) if np.isfinite(best)
                          else dict(a=[], chi2=[]))
        self.chi2_line.visible = bool(np.isfinite(best))
        if np.isfinite(best):
            self.chi2_line.location = best
        title = f"{ln.ion} {ln.lam0:.3f} Å"
        if a_ok:
            title += f"   {value_label} = {self.fmt(ln.a, ln.sigma_a, plus=True)}"
        if np.isfinite(best):
            title += f"   χ²ᵣ = {best:.2f}"
        self.fig.title.text = title
        self.crop(a_grid.tolist(), chi2.tolist(), self.fig)

    def clear(self):
        self.curve.data = dict(a=[], chi2=[])
        self.best.data = dict(a=[], chi2=[])
        self.fig.title.text = ""
        for layout in self._annotations:
            layout.visible = False
