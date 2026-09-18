"""The pair of Bokeh diagnostic panels every analysis result is read against.

An abundance, a correction and an ionisation balance all answer the same two questions of
their per-line values: is there a trend with excitation energy (a temperature diagnostic)
and one with reduced equivalent width (a microturbulence one). :class:`DiagnosticPanels`
owns that pair — two dark figures over one abundance axis, drawn from a single
:class:`~bokeh.models.ColumnDataSource` carrying an ``x_<axis>`` column per panel, so one
data push repaints both.

The points are linked to the per-line table in both directions without a server round trip:
a hovered point marks its row (:data:`HOVER_JS`) and a tapped point clicks it
(:data:`SELECT_JS`), through the ``idx`` column, which carries the index the row is
addressed by. The reverse direction (row → point) lives in the page script, which writes
the marker source published as ``<prefix>_hover``.

The panels carry no result of their own: the caller packs its lines into the source
(:meth:`DiagnosticPanels.pack_points`), fits the ranges to what is drawn
(:meth:`DiagnosticPanels.crop_panels`) and states the mean (:meth:`DiagnosticPanels.set_mean`).
A second series (the balance's ionised stage) is added with :meth:`DiagnosticPanels.add_series`.
"""
from __future__ import annotations

import numpy as np
from bokeh.models import (BoxAnnotation, ColumnDataSource, CustomJS, HoverTool, Range1d,
                          Span, Whisker)
from bokeh.plotting import figure

from tumult import abundance as ab

from ..ui.style import ACCENT, style_dark

# scatter x axes: key -> (axis label, per-line value fn)
X_AXES = {
    "evi": ("Excitation E_low [eV]", lambda ln: ln.e_low),
    "rew": ("Reduced EW  log(W/λ)", lambda ln: ab.reduced_ew(ln.ew, ln.lam0)),
}

# the two diagnostics every abundance fit is read against, shown as a pair of panels: the
# abundance against excitation energy (a temperature diagnostic) and against reduced EW (a
# microturbulence one).
PANEL_AXES = ("evi", "rew")

MEAN_COLOR = "#e05c5c"       # the mean line and its ± band

# what a hovered point says about itself
POINT_TOOLTIPS = [("λ₀", "@lam{0.000} Å"), ("ion", "@ion"), ("value", "@a{0.000}")]

# Hovering a point marks its table row: the point carries the index of the line it came from
# (the ``idx`` column), which is the same index the row is addressed by, so the two views are
# linked without a server round trip. The reverse direction (row → point) lives in the page
# script, which writes the marker source this module publishes as ``<prefix>_hover``.
HOVER_JS = """
document.querySelectorAll('.' + prefix + '-row.row-hover')
        .forEach((r) => r.classList.remove('row-hover'));
const ix = cb_data.index.indices;
const idx = source.data.idx;
if (!ix.length || !idx) return;
const row = document.querySelector(
    '.' + prefix + '-row[data-line="' + idx[ix[0]] + '"]');
if (row) row.classList.add('row-hover');
"""

# A scatter click selects through Bokeh's TapTool.  The source carries the original line
# index, so the callback can mark and reveal the corresponding table row, then click it to
# reuse the table's existing server-side jump action.
SELECT_JS = """
const selected = cb_obj.indices;
if (!selected.length || !source.data.idx) return;
const line = source.data.idx[selected[0]];
const table = document.getElementById(prefix + '-table');
if (!table) return;
const row = table.querySelector('.' + prefix + '-row[data-line="' + line + '"]');
if (!row) return;
// Reuse the row's ordinary click path for both the persistent table selection and the
// server-side line jump.
row.click();
const box = table.getBoundingClientRect();
const rect = row.getBoundingClientRect();
const head = table.querySelector('thead');
const top = box.top + (head ? head.getBoundingClientRect().height : 0);
if (rect.top < top || rect.bottom > box.bottom) {
    table.scrollTop += (rect.top + rect.height / 2) - (top + box.bottom) / 2;
}
"""


def empty_points():
    """The scatter columns with no line in them."""
    return ({f"x_{k}": [] for k in PANEL_AXES}
            | dict(a=[], lo=[], hi=[], lam=[], ion=[], idx=[]))


class DiagnosticPanels:
    """The two panels (:data:`PANEL_AXES`) of one analysis view, side by side.

    ``prefix`` namespaces the Bokeh model names the page script addresses
    (``<prefix>_points``, ``<prefix>_hover``, ``<prefix>_fig``/``<prefix>_fig2``);
    ``y_label`` is the initial abundance-axis label, which the view restates per result
    through :meth:`set_y_label`.
    """

    def __init__(self, prefix, y_label, source=None):
        self.prefix = prefix
        self.y_label = y_label
        self.figs, self.mean_lines, self.mean_bands = [], [], []
        self.src = source if source is not None else ColumnDataSource(empty_points())
        self.src.name = f"{prefix}_points"
        for i, key in enumerate(PANEL_AXES):
            fig = self._new_fig()
            fig.xaxis.axis_label = X_AXES[key][0]
            fig.name = f"{prefix}_fig" + ("" if i == 0 else str(i + 1))
        # the ring the page script moves onto the point of a hovered table row, on every
        # panel at once (it carries both panels' x, like the points themselves)
        self.hover_src = ColumnDataSource({f"x_{k}": [] for k in PANEL_AXES} | dict(a=[]))
        self.hover_src.name = f"{prefix}_hover"
        for key, fig in zip(PANEL_AXES, self.figs):
            fig.scatter(f"x_{key}", "a", source=self.hover_src, size=16, fill_color=None,
                        line_color="#fff", line_width=2)
        self.add_series(self.src, ACCENT)

    def _new_fig(self):
        fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                     tools="pan,wheel_zoom,box_zoom,tap,reset", active_scroll="wheel_zoom",
                     min_border_left=62, min_border_bottom=42,
                     x_range=Range1d(0.0, 1.0), y_range=Range1d(0.0, 1.0))
        style_dark(fig)   # match the other Bokeh panels' dark theme
        span = Span(location=0, dimension="width", line_color=MEAN_COLOR,
                    line_width=1.2, line_dash="dashed", visible=False)
        band = BoxAnnotation(fill_color=MEAN_COLOR, fill_alpha=0.08, visible=False)
        fig.add_layout(span)
        fig.add_layout(band)
        fig.yaxis.axis_label = self.y_label
        self.figs.append(fig)
        self.mean_lines.append(span)
        self.mean_bands.append(band)
        return fig

    def set_y_label(self, label):
        """Restate the abundance axis of both panels (what ``ln.a`` is for this result)."""
        for fig in self.figs:
            fig.yaxis.axis_label = label

    def add_series(self, source, color, xkeys=None):
        """Draw ``source``'s points on every panel, each linked to the per-line table."""
        source.selected.js_on_change(
            "indices", CustomJS(args=dict(source=source, prefix=self.prefix), code=SELECT_JS))
        for key, fig in zip(xkeys or [f"x_{k}" for k in PANEL_AXES], self.figs):
            r = fig.scatter(key, "a", source=source, size=8, fill_color=color,
                            line_color="#fff", line_width=1,
                            nonselection_fill_color=color, nonselection_fill_alpha=1.0,
                            nonselection_line_color="#fff", nonselection_line_alpha=1.0)
            fig.add_tools(HoverTool(
                renderers=[r], tooltips=POINT_TOOLTIPS,
                callback=CustomJS(args=dict(source=source, prefix=self.prefix),
                                  code=HOVER_JS)))

    def add_whiskers(self, source, color):
        """Error bars on both panels, from ``source``'s ``lo``/``hi`` columns."""
        for key, fig in zip(PANEL_AXES, self.figs):
            fig.add_layout(Whisker(base=f"x_{key}", upper="hi", lower="lo", source=source,
                                   line_color=color, line_width=1.2))

    @staticmethod
    def pack_points(indexed, value_of=lambda ln: ln.a, sigma_of=lambda ln: ln.sigma_a):
        """The scatter columns for ``indexed`` — ``(index into res.lines, line)`` pairs.

        The bars are part of the data, not decoration: they go into ``lo``/``hi`` so the
        ranges are fitted to what is drawn, including a point whose bar reaches furthest.
        The index rides along as ``idx``, which is what links a point to its table row.
        """
        lines = [ln for _i, ln in indexed]
        sig = [s if np.isfinite(s) else 0.0 for s in (float(sigma_of(ln)) for ln in lines)]
        vals = [float(value_of(ln)) for ln in lines]
        data = {f"x_{k}": [X_AXES[k][1](ln) for ln in lines] for k in PANEL_AXES}
        data.update(a=vals, lo=[v - s for v, s in zip(vals, sig)],
                    hi=[v + s for v, s in zip(vals, sig)],
                    lam=[ln.lam0 for ln in lines], ion=[ln.ion for ln in lines],
                    idx=[i for i, _ln in indexed])
        return data

    def crop(self, xs, ys, fig=None):
        """Fix a panel's ranges to the shown values with a 5% margin each side."""
        fig = fig if fig is not None else self.figs[0]
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        xs, ys = xs[np.isfinite(xs)], ys[np.isfinite(ys)]
        if xs.size == 0 or ys.size == 0:
            return
        for rng, lo, hi in ((fig.x_range, xs.min(), xs.max()),
                            (fig.y_range, ys.min(), ys.max())):
            pad = 0.05 * (hi - lo) or 0.05 * (abs(hi) or 1.0)
            start, end = float(lo - pad), float(hi + pad)
            rng.start, rng.end = start, end
            # Keep Reset useful after replacing Bokeh's auto-ranging DataRange1d with the
            # fixed range required for Whisker extents.
            if isinstance(rng, Range1d):
                rng.reset_start, rng.reset_end = start, end

    def crop_panels(self, series):
        """Fit every panel to its x data and the full vertical whisker span.

        Bokeh's automatic ranges do not account for ``Whisker`` annotations.  Derive the
        shared y extent here from the lower and upper error-bar columns so callers cannot
        accidentally crop to the point centres and clip a bar.
        """
        ys = [y for data in series for column in ("lo", "hi") for y in data[column]]
        for key, fig in zip(PANEL_AXES, self.figs):
            xs = [x for d in series for x in d[f"x_{key}"]]
            self.crop(xs, ys, fig)

    def set_mean(self, value, err=0.0):
        """Show the mean (and its ± band) on every panel, or hide it when there is none."""
        ok = value is not None and np.isfinite(value)
        err = err if np.isfinite(err) else 0.0
        for span, band in zip(self.mean_lines, self.mean_bands):
            span.visible = band.visible = bool(ok)
            if ok:
                span.location = float(value)
                band.top, band.bottom = float(value) + err, float(value) - err
