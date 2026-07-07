"""Bokeh view: an overview (nav) plot above a detailed (main) plot.

Holds any number of objects (loaded observations, tumult runs, synthetic spectra)
at once — observations share one colour and are distinguished by **line dash**
style (solid, dashed, dotted, …); runs get distinct colours. Windowing /
cropping happens in each spectrum's native units; only the decimated result is
converted to the chosen display units.
"""
from __future__ import annotations

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, DataRange1d, Range1d, RangeTool, Span, WheelZoomTool
from bokeh.plotting import figure

from tumult import Spectrum
from tumult.utils import pchip_mono

from debounce import Debouncer
from units import X_UNITS, Y_UNITS, label_for, to_value, to_range

# percentile of the flux used as the flat continuum when a spectrum is first loaded,
# so raw/continuum lands at ~0–1 without a single noise spike pinning the top
NORM_PERCENTILE = 98.0

# the nav overview autoscales its y-axis to the visible window, but the top is capped
# here so a noisy spike / emission peak can't squash the normalised continuum flat
NAV_Y_MAX = 1.05

# --- palette (mirrors STYLE_REFERENCE.md) ------------------------------------
BG = "#1e1e1e"
SURFACE = "#2a2a2a"
BORDER = "#3a3a3a"
GRID = "#2c2c2c"
MUTED = "#888888"
ACCENT = "#4a90e2"
# distinct from ACCENT so the abundance best-fit synthetic reads clearly against the (ACCENT)
# EW hatch it is drawn over
SYNTH_FIT = "#3ec46d"
LINE = "#cfd2d6"

# fixed line-colour palette; PALETTE[0] == LINE is the default (gray) object colour,
# PALETTE[1:] are the non-gray colours assigned to tumult runs. The first run colour is a
# cornflower blue kept distinct from ACCENT (#4a90e2) so a run never matches the UI accent.
PALETTE = ["#cfd2d6", "#7b8cde", "#50c878", "#e0a030", "#e05c5c", "#b084e8", "#34c9a0"]
RUN_COLORS = PALETTE[1:]

# objects are distinguished by dash; runs by colour (then dash once colours run out)
DASHES = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
DASH_CSS = {"solid": "solid", "dashed": "dashed", "dotted": "dotted",
            "dotdash": "dashed", "dashdot": "dotted"}


def _style(p):
    p.background_fill_color = BG
    p.border_fill_color = BG  # dark plot surround to match the darkened work surface
    p.outline_line_color = BORDER
    p.grid.grid_line_color = GRID
    p.axis.axis_line_color = BORDER
    p.axis.major_tick_line_color = BORDER
    p.axis.minor_tick_line_color = None
    p.axis.major_label_text_color = MUTED
    p.axis.axis_label_text_color = MUTED
    p.toolbar.logo = None


# public alias so other panels (e.g. the abundance scatter) get the same dark theme
style_dark = _style


class SpectrumObject:
    """One loaded spectrum and its renderers on the two plots."""

    def __init__(self, spectrum, name, dash, uid):
        self.spectrum = spectrum   # the displayed, continuum-normalised spectrum (~0–1)
        self.raw = spectrum        # the originally loaded spectrum, in its native units
        self.norm_knots = None     # (kx, ky) continuum-spline knots in native units, or None
        self.continuum = None      # last continuum evaluated over raw.xdata (for the overlay)
        self.already_normalised = False  # treat the spectrum as pre-normalised (continuum ≡ 1)
        self.name = name
        self.dash = dash
        self.uid = uid
        self.visible = True
        self.source = None  # (kind, path) where it was loaded from, for session save/restore
        self.meta = None  # provenance for the Info tab: FITS header cards or synth params
        self.filtered = None  # filtered Spectrum (reversible); used when canvas.show_filtered
        self.conv = {}  # convolution params applied to the main view (vshift/vdop/…)
        self.main_source = ColumnDataSource(dict(x=[], y=[]))
        self.nav_source = ColumnDataSource(dict(x=[], y=[]))
        self.line = self.markers = self.nav_line = None  # set when added to a canvas


class SpectrumCanvas:
    def __init__(self, spectrum=None, name="Object 1", init_width=25.0,
                 n_main=5000, n_nav=2500, marker_threshold=300, debounce_ms=50,
                 on_change=None):
        self.objects: list[SpectrumObject] = []
        self._next_uid = 0
        self.selected_uid = None
        self.show_filtered = False  # when True, plot each obj's filtered spectrum
        self.norm_edit_uid = None   # uid being continuum-edited (shown raw, others hidden)
        self.norm_edit_spectrum = None  # tool-supplied spectrum to draw for it (e.g. smoothed)
        self.n_main = n_main
        self.n_nav = n_nav
        self.marker_threshold = marker_threshold
        self._debounce = Debouncer(debounce_ms)  # coalesce range-driven redraws
        self.on_change = on_change  # called with the object list whenever it changes

        # x display unit starts at the first spectrum's native unit (Å when the canvas
        # is empty, until the first object is added); y is always the normalised
        # (dimensionless) axis since every spectrum is continuum-normalised
        self.xunit = spectrum.xunit if spectrum is not None else u.AA
        self.yunit = u.dimensionless_unscaled
        self.xlabel = label_for(self.xunit, X_UNITS, str(self.xunit))
        self.ylabel = label_for(self.yunit, Y_UNITS, str(self.yunit))

        # placeholder span for an empty canvas; real bounds are set on the first add
        if spectrum is not None:
            x, _ = spectrum.values
            wmin, wmax = float(x.min()), float(x.max())
        else:
            wmin, wmax = 0.0, init_width

        # --- main (detail) plot --------------------------------------------
        main = figure(
            sizing_mode="stretch_both", output_backend="webgl",
            tools="xpan,xwheel_zoom,xbox_zoom,reset", active_scroll="xwheel_zoom",
            x_range=Range1d(wmin, min(wmax, wmin + init_width), bounds=(wmin, wmax)),
            y_range=DataRange1d(range_padding=0.05, range_padding_units="percent",
                                only_visible=True),  # hidden objects drop out of the y-fit
            x_axis_label=self._x_axis_label(), y_axis_label=self._y_axis_label(),
            min_border_left=62, min_border_bottom=42,  # fixed frame insets (see .ew-controls)
        )
        _style(main)
        self.main = main

        # --- nav (overview) plot -------------------------------------------
        nav_wheel = WheelZoomTool(dimensions="width")
        nav = figure(height=96, sizing_mode="stretch_width", output_backend="webgl",
                     tools=[nav_wheel], active_scroll=nav_wheel, toolbar_location=None,
                     x_range=Range1d(wmin, wmax, bounds=(wmin, wmax)),
                     y_range=Range1d(0.0, NAV_Y_MAX),  # autofit each refresh, top capped
                     min_border_left=62)  # align the overview frame with the main frame
        nav.yaxis.visible = False
        nav.xaxis.axis_label = None
        _style(nav)
        range_tool = RangeTool(x_range=main.x_range, start_gesture="pan")
        range_tool.overlay.fill_color = ACCENT
        range_tool.overlay.fill_alpha = 0.12
        range_tool.overlay.line_color = ACCENT
        nav.add_tools(range_tool)
        self.nav = nav

        # faint copy of the original raw spectrum, drawn behind the edited one during
        # continuum normalisation so the smoothing/selection can be compared against it
        self.norm_raw_source = ColumnDataSource(dict(x=[], y=[]))
        self.norm_raw_line = main.line("x", "y", source=self.norm_raw_source,
                                       color=LINE, line_width=1, line_alpha=0.25,
                                       level="underlay", visible=False)

        # a synthetic line drawn over the observation while the abundance tool jumps to a
        # diagnostic line, so the fit quality is visible during manual line curation
        self.abund_overlay_source = ColumnDataSource(dict(x=[], y=[]))
        self.abund_overlay_line = main.line("x", "y", source=self.abund_overlay_source,
                                            color=SYNTH_FIT, line_width=1.5, line_dash="dashed",
                                            level="overlay", visible=False)

        # the fitted-profile integrated area for the jumped line, hatched exactly like the EW
        # tool so the EW measurement reviewed in step 1 reads the same on both surfaces
        self.abund_area_source = ColumnDataSource(dict(x=[], y1=[], y2=[]))
        self.main.varea(x="x", y1="y1", y2="y2", source=self.abund_area_source,
                        fill_color=ACCENT, fill_alpha=0.10,
                        hatch_pattern="/", hatch_color=ACCENT, hatch_alpha=0.5,
                        level="underlay")
        self.abund_area_line = self.main.line("x", "y1", source=self.abund_area_source,
                                              color=ACCENT, line_width=1.2, level="underlay",
                                              visible=False)

        # vertical bounds for the jumped line: the auto-crop fit window (amber) and the wider
        # synthesis window (red). Drawn between the discrete spectrum points so it is clear
        # which points fall in each window during line curation.
        self.abund_fit_spans = [
            Span(location=0, dimension="height", line_color="#e0a030", line_width=1.5,
                 line_dash="dashed", line_alpha=0.9, level="underlay", visible=False)
            for _ in range(2)]
        self.abund_synth_spans = [
            Span(location=0, dimension="height", line_color="#e05c5c", line_width=1.5,
                 line_dash="dashed", line_alpha=0.9, level="underlay", visible=False)
            for _ in range(2)]
        for _sp in (*self.abund_fit_spans, *self.abund_synth_spans):
            main.add_layout(_sp)

        main.x_range.on_change("start", self._on_main_range)
        main.x_range.on_change("end", self._on_main_range)
        nav.x_range.on_change("start", self._on_nav_range)
        nav.x_range.on_change("end", self._on_nav_range)

        # nav and main are embedded into separate template slots, so the area
        # above the main figure can swap between the nav overview and the line table
        self.main.name = "main_fig"
        self.nav.name = "nav_fig"

        if spectrum is not None:
            self.add_object(spectrum, name)

    # --- objects -----------------------------------------------------
    def next_run_style(self):
        """Pick a (color, dash) for a tumult run: the first unused non-gray palette
        colour, then iterate line styles once every colour is taken."""
        used = {(o.line.glyph.line_color, o.dash) for o in self.objects}
        used_colors = {c for c, _ in used}
        for color in RUN_COLORS:
            if color not in used_colors:
                return color, DASHES[0]
        for dash in DASHES:
            for color in RUN_COLORS:
                if (color, dash) not in used:
                    return color, dash
        n = len(self.objects)
        return RUN_COLORS[n % len(RUN_COLORS)], DASHES[n % len(DASHES)]

    def add_object(self, spectrum, name=None, color=None, dash=None, source=None, meta=None):
        first = not self.objects
        name = name or f"Object {len(self.objects) + 1}"
        color = color or LINE
        if first and spectrum.xunit != self.xunit:
            # the first object on an empty canvas sets the display x unit to its native one
            self.xunit = spectrum.xunit
            self.xlabel = label_for(self.xunit, X_UNITS, str(self.xunit))
            self.main.xaxis.axis_label = self._x_axis_label()
        if dash is None:
            # first line style not currently in use (cycle once every style is taken)
            used = {o.dash for o in self.objects}
            dash = next((d for d in DASHES if d not in used),
                        DASHES[len(self.objects) % len(DASHES)])
        uid = str(self._next_uid)
        self._next_uid += 1
        obj = SpectrumObject(spectrum, name, dash, uid)
        obj.source = source
        obj.meta = meta
        self.normalise_object(obj)  # raw -> ~0–1; flat continuum until a spline is set
        obj.line = self.main.line("x", "y", source=obj.main_source,
                                  color=color, line_width=1, line_dash=dash)
        obj.markers = self.main.scatter("x", "y", source=obj.main_source, size=5,
                                        fill_color=color, line_color=color,
                                        fill_alpha=0.55, line_alpha=0.55,
                                        level="underlay",  # always beneath the lines
                                        visible=False)
        obj.nav_line = self.nav.line("x", "y", source=obj.nav_source,
                                     color=color, line_width=1, line_dash=dash)
        self.objects.append(obj)
        self.selected_uid = obj.uid  # newly loaded becomes the selected/inspected one

        lo, hi = self.display_bounds()
        self.main.x_range.bounds = (lo, hi)
        self.nav.x_range.bounds = (lo, hi)
        self.nav.x_range.start, self.nav.x_range.end = lo, hi  # overview shows everything
        if first or (source and source[0] == "file"):
            # a newly loaded object, or the first object on a previously empty canvas:
            # focus the detail view on its wavelength range (the placeholder range that
            # an empty canvas carries has nothing worth preserving)
            self.main.x_range.start, self.main.x_range.end = self._display_range(spectrum)
        else:
            self.main.x_range.start = max(lo, self.main.x_range.start)
            self.main.x_range.end = min(hi, self.main.x_range.end)

        self._refresh_main()
        self._refresh_nav()
        if self.on_change:
            self.on_change(self)
        return obj

    def _find(self, uid):
        return next((o for o in self.objects if o.uid == uid), None)

    def selected(self):
        return self._find(self.selected_uid)

    def select_object(self, uid):
        if self._find(uid) is not None:
            self.selected_uid = uid
            if self.on_change:
                self.on_change(self)

    def toggle_object(self, uid):
        obj = self._find(uid)
        if obj is None:
            return
        obj.visible = not obj.visible
        obj.line.visible = obj.visible
        obj.nav_line.visible = obj.visible
        if obj.visible:
            self._refresh_object(obj)
        else:
            obj.markers.visible = False
            self._fit_nav_y()  # re-expand the nav y-fit now this object dropped out
        if self.on_change:
            self.on_change(self)

    def remove_object(self, uid):
        obj = self._find(uid)
        if obj is None:
            return
        self.main.renderers = [r for r in self.main.renderers
                               if r not in (obj.line, obj.markers)]
        self.nav.renderers = [r for r in self.nav.renderers if r is not obj.nav_line]
        self.objects.remove(obj)
        if self.selected_uid == uid:
            self.selected_uid = self.objects[0].uid if self.objects else None
        if self.objects:
            lo, hi = self.display_bounds()
            self.main.x_range.bounds = (lo, hi)
            self.nav.x_range.bounds = (lo, hi)
            # pull the ranges back inside the new bounds, else the tools (e.g.
            # nav wheel-zoom) lock up because the range sits outside its bounds
            self.nav.x_range.start, self.nav.x_range.end = lo, hi
            ms = max(lo, min(self.main.x_range.start, hi))
            me = min(hi, max(self.main.x_range.end, lo))
            self.main.x_range.start, self.main.x_range.end = (ms, me) if ms < me else (lo, hi)
            self._refresh_main()
            self._refresh_nav()
        if self.on_change:
            self.on_change(self)

    def display_bounds(self):
        ranges = [self._display_range(o.spectrum) for o in self.objects]
        return min(r[0] for r in ranges), max(r[1] for r in ranges)

    def selected_display_range(self):
        """Data extent (in the canvas xunit) of the active object, or the full
        extent across all objects when none is selected."""
        obj = self.selected()
        if obj is None:
            return self.display_bounds()
        return self._display_range(obj.spectrum)

    def _display_range(self, spectrum):
        x, _ = spectrum.values
        return to_range(float(x.min()), float(x.max()), spectrum.xunit, self.xunit)

    # --- axis labels -------------------------------------------------------
    def _x_axis_label(self):
        kind = "Energy" if self.xunit == u.eV else "Wavelength"
        return f"{kind} [{self.xlabel}]"

    def _y_axis_label(self):
        if self.yunit == u.dimensionless_unscaled:
            return "Normalized flux"
        return f"Flux [{self.ylabel}]"

    # --- unit switching ----------------------------------------------------
    def set_xunit(self, label):
        new, old = X_UNITS[label], self.xunit
        self.xunit, self.xlabel = new, label
        self.main.xaxis.axis_label = self._x_axis_label()
        for rng in (self.main.x_range, self.nav.x_range):
            self._reunit_range(rng, old, new)
        self._refresh_nav()
        self._refresh_main()

    def set_yunit(self, label):
        self.yunit, self.ylabel = Y_UNITS[label], label
        self.main.yaxis.axis_label = self._y_axis_label()
        self._refresh_nav()
        self._refresh_main()

    @staticmethod
    def _reunit_range(rng, old, new):
        if rng.bounds:
            rng.bounds = to_range(rng.bounds[0], rng.bounds[1], old, new)
        rng.start, rng.end = to_range(rng.start, rng.end, old, new)  # ascending across eV flip

    # --- debounced range handlers -----------------------------------------
    def _on_main_range(self, attr, old, new):
        self._debounce("main", self._refresh_main)

    def _on_nav_range(self, attr, old, new):
        self._debounce("nav", self._refresh_nav)

    def reset_nav_view(self):
        """Return the nav overview to the full spectral extent (the line navigator
        zooms it to a single line, so leaving the line view restores the overview)."""
        if not self.objects:
            return
        lo, hi = self.display_bounds()
        self.nav.x_range.start, self.nav.x_range.end = lo, hi

    # --- decimated redraws (window in native units, convert the result) ----
    def _window(self, spectrum, x_range, n, conv=None, convert_y=True):
        # convolution (crop kwargs) is cheap here: crop happens first, so only the
        # windowed slice is convolved — applied to both the main view and the nav
        a, b = to_range(x_range.start, x_range.end, self.xunit, spectrum.xunit)
        sub = spectrum.crop(xmin=a, xmax=b, **(conv or {})).decimate(n)
        try:
            # in continuum-edit mode the raw spectrum is shown in its native flux unit,
            # so only the x-axis is converted (convert_y=False)
            sub = sub.convert_units(xunit=self.xunit,
                                    yunit=self.yunit if convert_y else None)
        except Exception:
            pass  # incompatible units (e.g. flux vs normalized) — plot as-is
        return sub.values

    # --- continuum normalisation ------------------------------------------
    def _flat_continuum(self, raw):
        """Identity for already-normalised runs, else a flat robust-percentile level."""
        if raw is None:
            return None
        _, y = raw.values
        return np.full(y.shape, float(np.nanpercentile(y, NORM_PERCENTILE)))

    def normalise_object(self, obj):
        """(Re)build obj.spectrum = raw / continuum as a dimensionless (~0–1) spectrum.

        The continuum is the monotonic spline through obj.norm_knots when set, otherwise
        a flat level (identity for synthetic runs, which are already continuum-normalised).
        """
        raw = obj.raw
        x, y = raw.values
        if obj.already_normalised:
            cont = np.ones_like(y)   # assume the spectrum is already normalised (continuum ≡ 1)
        elif obj.norm_knots is not None:
            kx, ky = np.asarray(obj.norm_knots[0], float), np.asarray(obj.norm_knots[1], float)
            # clamp to the knot span so the spline holds flat past the ends rather than
            # diverging on the cubic extrapolation
            cont = pchip_mono(kx, ky)(np.clip(x, kx.min(), kx.max()))
        elif obj.source and obj.source[0] in ("run", "synth"):
            cont = np.ones_like(y)
        else:
            cont = self._flat_continuum(raw)
        cont = np.where(np.isfinite(cont) & (cont != 0.0), cont, np.nan)
        obj.continuum = cont
        obj.spectrum = Spectrum(raw.xdata, y / cont, yunit=u.dimensionless_unscaled,
                                xref=raw.xref)

    def renormalise(self, obj):
        """Recompute one obj's normalisation and redraw it (after its knots change)."""
        self.normalise_object(obj)
        if obj.uid != self.norm_edit_uid:  # edit view shows raw, so only redraw off-mode
            self._refresh_object(obj)

    # --- continuum-edit mode (show one spectrum raw, hide the rest) --------
    def _display(self, obj):
        """(spectrum, raw_mode) to plot for obj: the raw (or tool-smoothed) native flux
        while it is being continuum-edited."""
        if obj.uid == self.norm_edit_uid:
            return (self.norm_edit_spectrum or obj.raw), True
        return self.effective(obj), False

    def set_norm_view(self, spectrum):
        """Swap the spectrum shown for the edited obj (e.g. a low-pass-smoothed raw)."""
        self.norm_edit_spectrum = spectrum
        self._refresh_main()
        self._refresh_nav()

    def set_norm_edit(self, uid):
        """Enter/leave continuum-edit mode: show only `uid` in raw units (None to leave)."""
        self.norm_edit_uid = uid
        self.norm_edit_spectrum = None
        self.norm_raw_line.visible = uid is not None
        if uid is None:
            self.norm_raw_source.data = dict(x=[], y=[])
        for o in self.objects:
            shown = o.visible if uid is None else (o.uid == uid)
            o.line.visible = shown
            o.nav_line.visible = shown
            if not shown:
                o.markers.visible = False
        if uid is not None:
            raw = self._find(uid).raw
            self.main.yaxis.axis_label = (
                f"Flux [{label_for(raw.yunit, Y_UNITS, str(raw.yunit))}]")
        else:
            self.main.yaxis.axis_label = self._y_axis_label()
        self._refresh_main()
        self._refresh_nav()

    # --- abundance line overlay (synthetic best-fit drawn on the observation) ----
    def set_abund_overlay(self, spectrum, window):
        """Draw ``spectrum`` (a synthetic line) over ``window`` (Å) on the main figure.

        ``spectrum`` None or a non-finite window clears the overlay. The synthetic is cropped
        to the window and converted to the display units, so it lines up with the observation.
        """
        if spectrum is None or window is None or not np.isfinite(window[0]):
            return self.clear_abund_overlay()
        try:
            sub = spectrum.crop(xmin=window[0], xmax=window[1]).decimate(self.n_main)
            sub = sub.convert_units(xunit=self.xunit, yunit=self.yunit)
            x, y = sub.values
        except Exception:
            return self.clear_abund_overlay()
        self.abund_overlay_source.data = dict(x=np.asarray(x), y=np.asarray(y))
        self.abund_overlay_line.visible = True

    def clear_abund_overlay(self):
        self.abund_overlay_source.data = dict(x=[], y=[])
        self.abund_overlay_line.visible = False

    def set_abund_hatch(self, func, window):
        """Hatch the fitted-profile integrated area for the jumped diagnostic line.

        ``func`` is the fitted profile ``f(λ_Å) -> normalised flux`` (from ``measure_line``);
        the area between it and the continuum (=1) over ``window`` (Å) is filled and hatched
        with the same ``/`` pattern the EW tool uses. ``func`` None / bad window clears it.
        """
        if func is None or window is None or not np.isfinite(window[0]):
            return self.clear_abund_hatch()
        try:
            xg = np.linspace(float(window[0]), float(window[1]), 512)
            y1 = np.asarray(func(xg), dtype=float)
            xd = to_value(xg, u.AA, self.xunit)
        except Exception:
            return self.clear_abund_hatch()
        self.abund_area_source.data = dict(x=np.asarray(xd), y1=y1, y2=np.ones_like(xg))
        self.abund_area_line.visible = True

    def clear_abund_hatch(self):
        self.abund_area_source.data = dict(x=[], y1=[], y2=[])
        self.abund_area_line.visible = False

    @staticmethod
    def _snap_edge(xpts, value, upper):
        """Nudge a window edge into the gap between its last included and first excluded point.

        ``xpts`` are the sorted reference wavelengths (Å); ``value`` is the edge (Å). For a
        lower edge (``upper=False``) the line lands halfway between the point just below and
        the point just at/above; for an upper edge, between the point just at/below and just
        above. Falls back to ``value`` when there is no point on one side.
        """
        if xpts is None or len(xpts) == 0:
            return value
        if upper:
            inside, outside = xpts[xpts <= value], xpts[xpts > value]
        else:
            inside, outside = xpts[xpts >= value], xpts[xpts < value]
        if len(inside) == 0 or len(outside) == 0:
            return value
        near_in = inside[-1] if upper else inside[0]
        near_out = outside[0] if upper else outside[-1]
        return 0.5 * (float(near_in) + float(near_out))

    def set_abund_windows(self, fit_window, synth_window, xpts=None):
        """Draw the fit-window (amber) and synth-window (red) bounds for the jumped line.

        Each edge is snapped into the gap between the last included and first excluded
        reference point (``xpts``, Å) so the membership of each spectrum point is unambiguous,
        then converted to the display wavelength unit. A None / non-finite window hides its
        pair of lines.
        """
        xp = None
        if xpts is not None and len(xpts):
            xp = np.sort(np.asarray(xpts, dtype=float))
        for spans, win in ((self.abund_fit_spans, fit_window),
                           (self.abund_synth_spans, synth_window)):
            if win is None or not (np.isfinite(win[0]) and np.isfinite(win[1])):
                for sp in spans:
                    sp.visible = False
                continue
            lo = self._snap_edge(xp, float(win[0]), upper=False)
            hi = self._snap_edge(xp, float(win[1]), upper=True)
            a, b = to_range(lo, hi, u.AA, self.xunit)
            for sp, loc in zip(spans, (a, b)):
                sp.location = float(loc)
                sp.visible = True

    def clear_abund_window(self):
        for sp in (*self.abund_fit_spans, *self.abund_synth_spans):
            sp.visible = False

    # --- reversible filtering ---------------------------------------------
    def effective(self, obj):
        """The spectrum currently shown for obj (the filtered one when filtering is on)."""
        if self.show_filtered and obj.filtered is not None:
            return obj.filtered
        return obj.spectrum

    def apply_filter(self, mapping):
        """Show filtered spectra (uid -> filtered Spectrum); originals are kept."""
        for uid, spec in mapping.items():
            obj = self._find(uid)
            if obj is not None:
                obj.filtered = spec
        self.show_filtered = True
        self._refresh_main()
        self._refresh_nav()

    def clear_filter(self):
        """Revert to the unfiltered spectra."""
        self.show_filtered = False
        for obj in self.objects:
            obj.filtered = None
        self._refresh_main()
        self._refresh_nav()

    def set_filter_visible(self, on):
        """Show/hide the filtered spectra without discarding them (reversible toggle)."""
        self.show_filtered = bool(on)
        self._refresh_main()
        self._refresh_nav()

    def refresh_object(self, obj):
        """Redraw one object's main view (used after its convolution changes)."""
        spec, raw_mode = self._display(obj)
        conv = None if raw_mode else obj.conv  # raw view is the un-convolved continuum
        x, y = self._window(spec, self.main.x_range, self.n_main, conv,
                            convert_y=not raw_mode)
        obj.main_source.data = dict(x=x, y=y)
        shown = raw_mode or (obj.visible and self.norm_edit_uid is None)
        obj.markers.visible = shown and len(x) <= self.marker_threshold
        if raw_mode:  # keep the faint original raw in sync behind the smoothed/selected view
            rx, ry = self._window(obj.raw, self.main.x_range, self.n_main, convert_y=False)
            self.norm_raw_source.data = dict(x=rx, y=ry)

    def set_style(self, obj, color=None, width=None, dash=None, alpha=None):
        for glyph in (obj.line.glyph, obj.nav_line.glyph):
            if color is not None:
                glyph.line_color = color
            if width is not None:
                glyph.line_width = width
            if dash is not None:
                glyph.line_dash = dash
            if alpha is not None:
                glyph.line_alpha = alpha
        if dash is not None:
            obj.dash = dash
        if self.on_change:
            self.on_change(self)  # refresh the dash swatch in the list

    def _refresh_object(self, obj):
        self.refresh_object(obj)
        spec, raw_mode = self._display(obj)
        conv = None if raw_mode else obj.conv  # raw view is the un-convolved continuum
        xn, yn = self._window(spec, self.nav.x_range, self.n_nav, conv, convert_y=not raw_mode)
        obj.nav_source.data = dict(x=xn, y=yn)
        self._fit_nav_y()

    def _refresh_main(self):
        for obj in self.objects:
            if obj.uid != self.norm_edit_uid and not obj.visible:
                continue
            if self.norm_edit_uid is not None and obj.uid != self.norm_edit_uid:
                continue  # only the edited spectrum is drawn while editing its continuum
            self.refresh_object(obj)

    def _refresh_nav(self):
        for obj in self.objects:
            if obj.uid != self.norm_edit_uid and not obj.visible:
                continue
            if self.norm_edit_uid is not None and obj.uid != self.norm_edit_uid:
                continue
            spec, raw_mode = self._display(obj)
            conv = None if raw_mode else obj.conv  # raw view is the un-convolved continuum
            x, y = self._window(spec, self.nav.x_range, self.n_nav, conv, convert_y=not raw_mode)
            obj.nav_source.data = dict(x=x, y=y)
        self._fit_nav_y()

    def _fit_nav_y(self):
        """Autoscale the nav y-axis to the data drawn in the visible window, capping the
        top at NAV_Y_MAX for the normalised view (raw continuum-edit flux, which is far
        above 1, keeps a plain data fit so it isn't clipped)."""
        ys = [obj.nav_source.data["y"] for obj in self._nav_drawn()
              if len(obj.nav_source.data["y"])]
        if not ys:
            return
        lo = min(float(np.nanmin(y)) for y in ys)
        hi = max(float(np.nanmax(y)) for y in ys)
        if not (np.isfinite(lo) and np.isfinite(hi)):
            return
        pad = 0.05 * ((hi - lo) or abs(hi) or 1.0)
        raw_mode = self.norm_edit_uid is not None  # only the raw-flux object is shown then
        top = hi + pad if raw_mode else min(hi + pad, NAV_Y_MAX)
        self.nav.y_range.start, self.nav.y_range.end = lo - pad, top

    def _nav_drawn(self):
        """Objects currently drawn in the nav (visible, or the lone edited one)."""
        for obj in self.objects:
            if self.norm_edit_uid is not None:
                if obj.uid == self.norm_edit_uid:
                    yield obj
            elif obj.visible:
                yield obj
