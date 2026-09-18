"""Spectrum filter tool with two selectable methods.

Clicking the toolbar Filter button opens the Filter inspector panel, whose method
selector switches between:

* **Optimal (Brault)** — the main figure is swapped for a "feautrier" view: the
  Fourier power spectrum (normalized frequency vs power, log y) of every loaded
  spectrum. A draggable crosshair sets the fit's ``fmax`` (its x) and the constant
  noise floor (its y, in log-power). The fit runs on the *selected* spectrum.
* **Low-pass** — a simple Gaussian low-pass in ln λ space (``Spectrum.lowpass``),
  parametrised by a single cutoff resolution ``R = λ/Δλ``. No fit / feautrier view.

Either way the resulting filter is applied to *all* spectra (reversibly — originals
are kept and the filter can be cleared / toggled).
"""
from __future__ import annotations

import logging
from enum import StrEnum

import numpy as np
from bokeh.models import (ColumnDataSource, CustomJS, PointDrawTool, Range1d, Span,
                          WheelZoomTool)
from bokeh.plotting import figure

from tumult.utils import (R_MAX, R_MIN, apply_brault_filter, apply_brault_popt,
                          brault_filter_H, brault_power_spectrum, brault_prep as _prep,
                          clamp_resolution)

from ..ui.style import ACCENT, LINE, MUTED, style_dark
from ..ui.controls import Field, Form
from ..ui.debounce import Debouncer
from ..ui.htmlbridge import pulse, toggle
from ..ui.latestwork import LatestWork
from tumult.spectrum_io import Spectrum
from ..ui.tool import Tool, suppress_callbacks

N_DISPLAY = 4000   # decimate power spectra for plotting
RED = "#e05c5c"
GREEN = "#50c878"  # optimal filter H

# The physical Fourier frequency (cycles per unit ln λ) IS the spectral resolution:
# the Nyquist frequency 1/(2·Δlnλ) is the max R a 2-pixel feature can have. So the
# axis value is R directly — no conversion, just a relabel capped at R_MAX. The signal
# preparation and those bounds are :mod:`tumult.utils`' (``brault_prep``,
# ``clamp_resolution``); the axis they label is this panel's.


class Method(StrEnum):
    """Filter method: the value is the internal key, ``label`` what the panel select shows."""
    BRAULT = "brault"
    LOWPASS = "lowpass"

    @property
    def label(self):
        return {Method.BRAULT: "Optimal (Brault)", Method.LOWPASS: "Low-pass"}[self]

    @classmethod
    def from_label(cls, label):
        """The method whose select label is ``label`` (defaults to Brault)."""
        return next((m for m in cls if m.label == label), cls.BRAULT)


class FilterTool(Tool):
    name = "filter"

    def __init__(self, ctx, debounce_ms=200):
        Tool.__init__(self, ctx)
        self.own_panel("filter-panel-content")
        self.own_mode_bridges((".sidebar-right", "filter-active"), ("#filter-btn", "active"))
        self._view_bridge = ctx.bridges.cls("#main", "filter-mode")   # feautrier view swap
        self._has_bridge = ctx.bridges.cls("#main", "has-filter")     # reveal the switch
        self._check_bridge = ctx.bridges.check("filter-switch-input")  # switch on/off state
        self.method = Method.LOWPASS    # resolution cut (LOWPASS) or feautrier fit (BRAULT)
        self.fmax = 1.0                 # physical frequency (cycles / unit ln λ); set on activate
        self.floor = -20.0              # log-power
        self._configured = False        # have we picked the default crosshair yet?
        self._suppress = False
        self._debounce = Debouncer(debounce_ms)  # coalesce crosshair-drag refits
        immediate = getattr(ctx.doc, "session_context", None) is None
        self._filter_work = LatestWork(lambda *args: ctx.submit(*args), "filter spectrum",
                                       lambda exc: self.notify(f"Filtering failed: {exc}"),
                                       immediate=immediate)
        self._power_work = LatestWork(lambda *args: ctx.submit(*args), "build filter view",
                                      lambda exc: self.notify(f"Filter view failed: {exc}"),
                                      immediate=immediate)

        # inspector panel: method selector + (low-pass only) cutoff resolution
        self.form = Form("filter", [
            Field("method", "Method", kind="select", value=Method.LOWPASS.label,
                  options=[m.label for m in Method]),
            Field("res", "Cutoff resolution R = λ/Δλ", kind="number",
                  value=20000.0, min=0),
        ], on_change=self._on_field)

        fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                     x_range=Range1d(R_MIN, R_MAX), y_range=Range1d(1e-12, 10.0),
                     x_axis_type="log", y_axis_type="log",
                     x_axis_label="resolution  R = λ / Δλ",
                     y_axis_label="power",
                     tools="xpan,xwheel_zoom,reset", active_scroll="xwheel_zoom",
                     min_border_left=62, min_border_bottom=42)
        style_dark(fig)
        # let the wheel widen the view past a range bound instead of sticking on the cursor
        for tool in fig.select(WheelZoomTool):
            tool.maintain_focus = False
        fig.name = "feautrier_fig"
        self.fig = fig

        self._sources = {}   # uid -> power-spectrum CDS
        self._lines = {}     # uid -> renderer

        # overlays drawn above the spectral data (level="annotation"), both in power:
        # the fitted brault_filter model, and model*H (the filtered/kept signal)
        self.model_source = ColumnDataSource(dict(x=[], y=[]))
        self.fig.line("x", "y", source=self.model_source, color=ACCENT, line_width=3,
                      level="annotation")
        self.h_source = ColumnDataSource(dict(x=[], y=[]))
        self.fig.line("x", "y", source=self.h_source, color=GREEN, line_width=3,
                      level="annotation")

        # draggable crosshair: x -> fmax, y -> noise floor (power); dashed guides
        self.cross = ColumnDataSource(dict(x=[self.fmax], y=[float(np.exp(self.floor))]))
        pt = self.fig.scatter("x", "y", source=self.cross, size=12, fill_color=RED,
                              line_color="#fff", line_width=1.5, level="overlay")
        draw = PointDrawTool(renderers=[pt], add=False, num_objects=1)
        self.fig.add_tools(draw)
        self.fig.toolbar.active_drag = draw
        self.vline = Span(location=self.fmax, dimension="height",
                          line_color=RED, line_dash="dashed", line_alpha=0.6)
        self.hline = Span(location=float(np.exp(self.floor)), dimension="width",
                          line_color=RED, line_dash="dashed", line_alpha=0.6)
        self.fig.add_layout(self.vline)
        self.fig.add_layout(self.hline)
        # guides track the dot live (client-side) while the server refit debounces
        self.cross.js_on_change("data", CustomJS(
            args=dict(v=self.vline, h=self.hline), code="""
            const d = cb_obj.data;
            if (d['x'].length) { v.location = d['x'][0]; h.location = d['y'][0]; }
        """))
        self.cross.on_change("data", self._on_drag)

    @property
    def models(self):
        return [*self.form.models]

    # --- activation -------------------------------------------------------
    def activate(self):
        sel = self.canvas.spectra.selected()
        if sel is None:
            self.notify("Load and select a spectrum before using the filter.")
            return
        self.active = True
        self.push_panel()
        self.push_mode("on")
        self._enter_method()

    def _enter_method(self):
        """Set up the view + apply the current method (called on activate / switch)."""
        if self.method == Method.LOWPASS:
            pulse(self._view_bridge, "off")  # low-pass uses the normal canvas, not feautrier
            self._apply_lowpass()
        else:
            self._enter_brault()

    def _switch_method(self, method):
        """Change the active method: drop the old filter, set up the new one."""
        if method == self.method:
            return
        self.method = method
        self.canvas.clear_filter()       # discard the previous method's result
        if method == Method.BRAULT:
            self._configured = False     # re-pick the crosshair for the fresh fit
        self._enter_method()
        self.push_panel()

    def _enter_brault(self):
        sel = self.canvas.spectra.selected()
        if sel is None:
            return
        spectra = [(o.uid, o.spectrum, o.uid == self.canvas.spectra.selected_uid)
                   for o in self.canvas.spectra]
        selected_uid = sel.uid

        def compute(token):
            powers = []
            for uid, spectrum, selected in spectra:
                token.check()
                _, y0, w, d = _prep(spectrum)
                s, npow = brault_power_spectrum(y0 * w, d)
                m = (s >= R_MIN) & (s <= R_MAX)
                powers.append((uid, selected, s[m], npow[m], d))
            return powers

        def apply(powers):
            if (not self.active or self.method != Method.BRAULT
                    or self.canvas.spectra.selected_uid != selected_uid
                    or [o.uid for o in self.canvas.spectra] != [row[0] for row in powers]):
                return
            self._show_power_spectra(powers)
            selected = next(row for row in powers if row[0] == selected_uid)
            self._finish_enter_brault(selected, powers)

        self._power_work.request(compute, apply)

    def _finish_enter_brault(self, selected, powers):
        """Finish entering the Brault view after its immutable power data arrive."""
        _, _, s, npow, d = selected
        logp = np.log(np.clip(npow, 1e-300, None))
        nyq = float(s.max())
        if not self._configured:
            # pick a default crosshair only once, so re-opening the panel keeps the
            # marker where the user (or the last fit) left it
            peak = float(logp.max())
            hi = logp[s > 0.5 * nyq]
            hi = hi[hi > peak - 690.0]          # ignore exact-zero (clipped) bins
            floor = float(np.median(hi)) if hi.size else peak - 6.0
            self.floor = min(floor, peak - 2.0)  # keep below the signal peak
            # default fmax: first frequency where the power has fallen to within
            # ~1 dex of the noise floor — the signal/noise crossover
            cross = (s > 0.02 * nyq) & (logp < self.floor + 1.0)
            self.fmax = float(np.clip(s[cross][0] if cross.any() else 0.3 * nyq, R_MIN, R_MAX))
            self._configured = True
        # log resolution axis, limited to [R_MIN, R_MAX]
        xmax = max(0.5 / abs(row[4]) for row in powers)
        self.fig.x_range.start, self.fig.x_range.end = R_MIN, min(nyq, R_MAX)
        self.fig.x_range.bounds = (R_MIN, min(max(xmax, nyq), R_MAX))
        # fixed y-range: with auto-range the marker drag rescales the axis (jumpy).
        # Anchor ~2 decades below the noise floor so the rolloff stays visible.
        allp = np.concatenate([np.asarray(src.data["y"]) for src in self._sources.values()])
        ymax = float(np.nanmax(allp)) * 2.0
        ymin = min(float(np.exp(self.floor)) * 1e-2, ymax * 1e-12)
        self.fig.y_range.start, self.fig.y_range.end = ymin, ymax
        self._set_cross(self.fmax, float(np.exp(self.floor)))
        self._refit_apply()
        pulse(self._view_bridge, "on")

    def deactivate(self):
        """Leave the filter mode; the filter stays applied (clear() removes it)."""
        self.active = False
        self._power_work.invalidate()
        self._filter_work.invalidate()
        pulse(self._view_bridge, "off")
        self.push_mode("off")

    def clear(self):
        self.canvas.clear_filter()
        self._configured = False
        self._emit_state()
        self.deactivate()

    def set_visible(self, on):
        """Toggle the filtered view in the canvas on/off, keeping the fitted filter."""
        self.canvas.set_filter_visible(on)
        self._emit_state()

    def _emit_state(self):
        """Reflect filter existence/visibility in the header switch, and tell the other
        tools (EW, abundance wizard) that the effective spectra changed."""
        exists = any(o.filtered is not None for o in self.canvas.spectra)
        visible = self.canvas.show_filtered
        toggle(self._has_bridge, exists)
        toggle(self._check_bridge, visible)
        self.ctx.emit("filter-state", exists=exists, visible=visible)
        self.ctx.workspace_changed()

    def update(self):
        """From main.update_all: keep the view in sync when the workspace changes."""
        if not self.active:
            return
        if self.method == Method.LOWPASS:
            self._apply_lowpass()
        else:
            self._enter_brault()

    # --- inspector panel --------------------------------------------------
    def _on_field(self, key, value):
        if not self.active:
            return
        if key == "method":
            self._switch_method(Method.from_label(value))
        elif key == "res" and self.method == Method.LOWPASS:
            self._debounce("lowpass", self._apply_lowpass)

    def panel_html(self):
        if not self.active:
            return ""
        f = self.form.fields
        f["res"].disabled = self.method != Method.LOWPASS
        parts = [f["method"].render("filter")]
        if self.method == Method.LOWPASS:
            parts.append(f["res"].render("filter"))
            hint = ("A Gaussian low-pass in ln λ space: suppresses features sharper than "
                    "the cutoff resolution. Simple and ringing-free — no fit required.")
        else:
            hint = ("Optimal Brault filter: drag the crosshair in the feautrier view "
                    "(x = max resolution, y = noise floor) to fit the signal/noise rolloff.")
        parts.append('<button class="text-btn filter-clear-btn">Clear filter</button>')
        parts.append(f'<div class="hint">{hint}</div>')
        return "".join(parts)

    # --- low-pass method --------------------------------------------------
    def _apply_lowpass(self):
        """Apply the simple resolution-cut low-pass to all spectra (interactive)."""
        if not self.active or self.method != Method.LOWPASS:
            return
        res = float(self.form.value("res"))
        if res <= 0:
            self.canvas.clear_filter()
            self._emit_state()
            return
        self._lowpass_filtered(res)

    def _lowpass_filtered(self, res):
        spectra = [(obj.uid, obj.spectrum) for obj in self.canvas.spectra]
        uids = [uid for uid, _ in spectra]

        def compute(token):
            out = {}
            for uid, spectrum in spectra:
                token.check()
                out[uid] = spectrum.lowpass(res)
            return out

        def apply(filtered):
            if [obj.uid for obj in self.canvas.spectra] != uids:
                return
            self.canvas.apply_filter(filtered)
            self._emit_state()

        self._filter_work.request(compute, apply)

    # --- power-spectrum display -------------------------------------------
    def _show_power_spectra(self, powers):
        """Install already-computed power spectra on the document thread."""
        self.fig.renderers = [r for r in self.fig.renderers
                              if r not in self._lines.values()]
        self._sources.clear()
        self._lines.clear()
        for uid, sel, s, npow, _ in powers:
            step = max(1, len(s) // N_DISPLAY)
            src = ColumnDataSource(dict(x=s[::step], y=npow[::step]))
            r = self.fig.line("x", "y", source=src, color=LINE if sel else MUTED,
                              line_width=1.5 if sel else 1, line_alpha=1.0 if sel else 0.5)
            self._sources[uid] = src
            self._lines[uid] = r

    # --- crosshair drag -> debounced refit --------------------------------
    @suppress_callbacks
    def _set_cross(self, x, y):
        self.cross.data = dict(x=[x], y=[y])
        self.vline.location = x
        self.hline.location = y

    def _on_drag(self, attr, old, new):
        if self._suppress or not self.active:
            return
        self._debounce("refit", self._refit_apply)

    # --- fit (selected) + apply (all) -------------------------------------
    def _refit_apply(self):
        """Fit + apply from the current crosshair position (interactive)."""
        if not self.active:
            return
        xs, ys = self.cross.data["x"], self.cross.data["y"]
        sel = self.canvas.spectra.selected()
        if not xs or not ys or sel is None:
            return
        fmax = clamp_resolution(xs[0], _prep(sel.spectrum)[3])
        floor = float(np.log(max(float(ys[0]), 1e-300)))
        self._apply(fmax, floor)

    def apply_saved(self, fmax, floor):
        """Re-apply a saved Brault filter (from a session) without opening the view."""
        if not self.canvas.spectra:
            return
        self.method = Method.BRAULT
        self.form.set_values(method=Method.BRAULT.label)
        self._configured = True
        self._set_cross(float(fmax), float(np.exp(floor)))
        try:
            self._apply(float(fmax), float(floor))
        except Exception:
            logging.exception("Restoring saved filter failed")

    def apply_saved_lowpass(self, res):
        """Re-apply a saved low-pass filter (from a session) without opening the view."""
        if not self.canvas.spectra:
            return
        self.method = Method.LOWPASS
        self.form.set_values(method=Method.LOWPASS.label, res=float(res))
        try:
            self._lowpass_filtered(float(res))
        except Exception:
            logging.exception("Restoring saved low-pass filter failed")

    def _apply(self, fmax, floor):
        self.fmax, self.floor = fmax, floor
        self.vline.location = fmax
        self.hline.location = float(np.exp(floor))
        sel = self.canvas.spectra.selected()
        if sel is None:
            return
        selected_uid = sel.uid
        spectra = [(obj.uid, obj.spectrum) for obj in self.canvas.spectra]
        uids = [uid for uid, _ in spectra]

        def compute(token):
            selected = next(spectrum for uid, spectrum in spectra if uid == selected_uid)
            _, y0_sel, w_sel, d_sel = _prep(selected)
            token.check()
            _, s, _, model, popt, _ = apply_brault_filter(
                y0_sel * w_sel, d_sel, fmax=fmax, noise_floor=floor)
            m = (s >= R_MIN) & (s <= R_MAX)
            sm, modm = s[m], model[m]
            Hm = brault_filter_H(sm, popt)
            filtered = {}
            for uid, spectrum in spectra:
                token.check()
                xx, y0, w, d = _prep(spectrum)
                y_smooth = apply_brault_popt(y0 * w, d, popt)
                fy = y_smooth + (1.0 - w) * y0 + 1.0
                filtered[uid] = Spectrum(xx, fy, xunit=spectrum.xunit,
                                         yunit=spectrum.yunit,
                                         xref=getattr(spectrum, "xref", None))
            return sm, modm, Hm, filtered

        def apply(result):
            if [obj.uid for obj in self.canvas.spectra] != uids:
                return
            sm, modm, hm, filtered = result
            step = max(1, len(sm) // N_DISPLAY)
            self.model_source.data = dict(x=sm[::step],
                                          y=np.clip(modm[::step], 1e-300, None))
            self.h_source.data = dict(x=sm[::step],
                                      y=np.clip((modm * hm)[::step], 1e-300, None))
            self.canvas.apply_filter(filtered)
            self._emit_state()

        self._filter_work.request(compute, apply)
