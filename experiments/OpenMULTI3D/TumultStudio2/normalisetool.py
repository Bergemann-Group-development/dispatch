"""Interactive spline continuum normalisation.

Clicking the toolbar **Normalise** button enters an edit mode: the canvas shows only
the selected spectrum in its raw (un-normalised) units, with a draggable continuum
spline drawn on top. The inspector swaps to a Normalise panel that can *initialise*
the spline (auto-place knots at the flux peaks) and tune the peak spacing. Knots are
added / dragged / deleted directly on the canvas (a ``PointDrawTool``); the spectrum
is divided by the monotonic spline through them, so leaving the mode shows the
continuum-normalised spectrum (continuum at ~1).

The mode is an edit session (:class:`tool.Tool`): entering it snapshots the spectrum's
knots, already-normalised flag and normalised spectrum. **Done** keeps the edits;
**Cancel** and Esc restore the snapshot. **Clear spline** removes the knots within the
session, returning the spectrum to the automatic flat scaling, which Cancel undoes too.
"""
from __future__ import annotations

import logging
import html

import astropy.units as u
import numpy as np
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, PointDrawTool
from scipy.signal import argrelextrema

from tumult import Spectrum, continuum
from tumult.utils import pchip_mono

from canvas import ACCENT, NORM_PERCENTILE
from controls import Field, Form
from debounce import Debouncer
from latestwork import LatestWork
from tool import Tool, suppress_callbacks
from units import to_value

RED = "#e05c5c"
N_SPLINE = 400   # samples of the drawn continuum curve


class NormaliseTool(Tool):
    name = "normalise"
    cancellable = True

    def __init__(self, ctx, debounce_ms=200):
        Tool.__init__(self, ctx)
        self.main = self.canvas.main
        self.own_panel("norm-panel-content")
        self.own_mode_bridges((".sidebar-right", "norm-active"), ("#normalise-btn", "active"))
        self._uid = None            # the object being edited
        self._suppress = False      # guard programmatic knot writes against the refit
        self._saved_range = None    # main x-range to restore when the tool closes
        self._selection_request = 0
        self._fit_spectrum = None   # rendered spectrum the spline is fit to (smoothed+selected)
        self._debounce = Debouncer(debounce_ms)
        self._normalise_work = LatestWork(
            lambda *args: ctx.submit(*args), "normalise continuum",
            lambda exc: self.notify(f"Continuum normalisation failed: {exc}"),
            immediate=getattr(ctx.doc, "session_context", None) is None)

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
        ], on_change=self._on_field, groups=[("Continuum detection settings", ("res", "win", "pct", "peak"))])

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
        self.begin_edit()
        self._seed_knots(obj)
        self.form.fields["prenorm"].value = obj.already_normalised
        self.canvas.set_norm_edit(obj.uid)   # show this spectrum raw, hide the others
        self._zoom_full()                    # fit the whole spectrum for continuum work
        self._sync_prenorm(obj)              # render + show/hide the spline per prenorm
        self.push_panel()
        self.push_mode("on")

    def deactivate(self):
        if not self.active:
            return
        if self._edit_baseline is not None and self._debounce.pending("refit"):
            self._refit()               # Done: apply a knot drag the debounce still holds
        self.end_edit()
        self.active = False
        self._uid = None
        self.spline_line.visible = False
        self.knot_glyph.visible = False
        self.canvas.release_edit_tool(self.draw)
        self.canvas.set_norm_edit(None)      # back to the normalised view
        self._restore_range()
        self.push_mode("off")

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

    def capture_edit(self):
        obj = self._obs()
        return dict(norm_knots=obj.norm_knots, already_normalised=obj.already_normalised,
                    continuum=obj.continuum, spectrum=obj.spectrum)

    def restore_edit(self, baseline):
        """Put the entry state back. Every edit replaces these attributes rather than
        mutating them, so the snapshot holds the very arrays the spectrum had."""
        obj = self._obs()
        if obj is None:
            return
        self._normalise_work.invalidate()
        for key, value in baseline.items():
            setattr(obj, key, value)
        self.form.fields["prenorm"].value = obj.already_normalised

    def clear(self):
        """Drop the spline: the spectrum returns to the automatic flat normalisation
        (not to the state on entry, which is what Cancel restores)."""
        obj = self._obs()
        if obj is None:
            return
        obj.norm_knots = None
        self._set_knots([], [])
        self._recompute(obj)                # back to the flat percentile continuum
        self._render_selection()
        self._redraw_spline(obj)                 # show the flat initial-continuum guide
        self.push_panel()

    def _sync_prenorm(self, obj):
        """Render the view and show/hide the spline affordances for obj.already_normalised."""
        on = obj.already_normalised
        self.spline_line.visible = not on
        self.knot_glyph.visible = not on
        if on:
            self.canvas.release_edit_tool(self.draw)
        else:
            self.canvas.claim_edit_tool(self.draw)
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
        self._recompute(obj)
        self._sync_prenorm(obj)
        self.push_panel()

    @staticmethod
    def _build_fit_spectrum(raw, settings):
        """The spectrum the spline is fit to and rendered — see
        :func:`tumult.continuum.reject_outliers`."""
        return continuum.reject_outliers(raw, resolution=settings["res"],
                                         window=settings["win"],
                                         percentile=settings["pct"])

    def _render_selection(self):
        self._prepare_selection(place_knots=False)

    def _prepare_selection(self, place_knots):
        """Compute from snapshots; only the latest request may update the canvas."""
        obj = self._obs()
        if obj is None or obj.already_normalised:
            return
        raw = obj.raw
        settings = {k: self.form.value(k) for k in ("res", "win", "pct", "peak")}
        self._selection_request += 1
        request = self._selection_request
        # Knot edits during a background fit take precedence over its result.
        knots = {k: list(v) for k, v in self.knots.data.items()}

        def compute(token):
            token.check()
            spec = self._build_fit_spectrum(raw, settings)
            ex = None
            if place_knots:
                x, y = spec.values
                if len(x) < 2:
                    raise ValueError("Not enough points to place continuum knots")
                dx = float(np.mean(np.diff(x)))
                nstep = max(1, int(abs(settings["peak"]) / abs(dx))) if dx else 1
                ex = argrelextrema(y, np.greater_equal, order=nstep)[0]
                ex = np.unique(np.concatenate(([0], ex, [len(x) - 1])))
            token.check()
            return spec, ex

        def apply(result):
            if (request != self._selection_request or not self.active
                    or self._obs() is not obj or obj.already_normalised
                    or any(self.form.value(k) != v for k, v in settings.items())):
                return
            spec, ex = result
            self._fit_spectrum = spec
            self.canvas.set_norm_view(spec)
            if ex is not None and knots == {k: list(v) for k, v in self.knots.data.items()}:
                x, y = spec.values
                kx = to_value(x[ex], spec.xunit, self.canvas.xunit)
                self._set_knots(kx, y[ex])
                self._recompute(obj)
                self._force_repaint()
                if ex.size <= 2:
                    self.notify("Few continuum peaks found — widen the peak spacing or lower the "
                                "resolution cut.", level="info")

        self.ctx.submit("normalise-peaks" if place_knots else "normalise-selection", compute, apply)

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
        """Place continuum knots without blocking the document callback."""
        self._prepare_selection(place_knots=True)

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
        fields_html = self.form.html()
        return (
            f'<div class="hint">Editing continuum: {html.escape(self._obs().name)}. Other spectra are temporarily hidden.</div>'
            f'{fields_html}'
            f'<button class="primary-btn norm-init"{dis}>Place continuum automatically</button>'
            f'<button class="text-btn norm-clear"{dis} title="Remove the knots: automatic flat '
            f'scaling, until Done (Cancel restores the entry state)">Clear spline</button>'
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
            self.push_panel()

    def _recompute(self, obj):
        """Rebuild the full-resolution normalisation away from the document thread."""
        d = self.knots.data
        kx_disp = np.asarray(d.get("x", []), float)
        ky = np.asarray(d.get("y", []), float)
        kx_nat = to_value(kx_disp, self.canvas.xunit, obj.raw.xunit)
        order = np.argsort(kx_nat)
        kx_nat, ky = kx_nat[order], ky[order]
        uniq = np.diff(kx_nat, prepend=-np.inf) > 0   # pchip needs increasing x
        kx_nat, ky = kx_nat[uniq], ky[uniq]
        if kx_nat.size < 2:
            knots = None
        else:
            knots = (kx_nat, ky)
        raw = obj.raw
        source = obj.source
        already_normalised = obj.already_normalised
        uid = obj.uid
        self._redraw_knots(raw, knots, obj.continuum)

        def compute(token):
            x, y = raw.values
            token.check()
            if already_normalised:
                cont = np.ones_like(y)
            elif knots is not None:
                kx, kvals = continuum.merge_knots(
                    np.asarray(x, float), np.asarray(y, float), base_knots=knots,
                    percentile=NORM_PERCENTILE)
                cont = continuum.evaluate(kx, kvals, x)
            elif source and source[0] in ("run", "synth"):
                cont = np.ones_like(y)
            else:
                cont = continuum.flat_continuum(y, NORM_PERCENTILE)
            cont = continuum.sanitise(cont)
            spectrum = Spectrum(raw.xdata, y / cont,
                                yunit=u.dimensionless_unscaled, xref=raw.xref)
            return knots, cont, spectrum

        def apply(result):
            current = self.canvas._find(uid)
            if current is not obj or current.raw is not raw:
                return
            current.norm_knots, current.continuum, current.spectrum = result
            if current.uid != self.canvas.norm_edit_uid:
                self.canvas._refresh_object(current)
            self._redraw_spline(current)
            self.ctx.workspace_changed()

        self._normalise_work.request(compute, apply)

    def _redraw_knots(self, raw, knots, flat):
        """Cheap spline preview while the full spectrum is being normalised."""
        x, _ = raw.values
        if knots is None:
            if flat is None or len(x) < 2:
                self.spline_source.data = dict(x=[], y=[])
                return
            xg = np.array([x[0], x[-1]])
            yg = np.array([flat[0], flat[-1]])
        else:
            kx_nat, ky = knots
            xg = np.linspace(kx_nat.min(), kx_nat.max(), N_SPLINE)
            yg = pchip_mono(kx_nat, ky)(xg)
        self.spline_source.data = dict(
            x=to_value(xg, raw.xunit, self.canvas.xunit), y=yg)

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

    @suppress_callbacks
    def _set_knots(self, kx, ky):
        self.knots.data = dict(x=list(map(float, kx)), y=list(map(float, ky)))
