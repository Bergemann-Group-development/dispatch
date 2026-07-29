"""The shared machinery of the two 3-step abundance wizards.

:class:`abundancetool.AbundanceTool` (element abundances vs an observation) and
:class:`abundcorrtool.AbundCorrTool` (NLTE-LTE corrections vs a synthetic 1D LTE
reference) grew as siblings with near-verbatim copies of the step switcher, the
results-list / selection management, the action-channel dispatch, the per-line table
scaffolding, the dark scatter figure and the jump-to-line logic. That common machinery
lives here; the subclasses keep what is genuinely different - their forms, their run
stages (what gets synthesised and inverted) and their table columns.

Subclasses set the class attributes below (``prefix`` namespaces the DOM/bridge ids,
``kind`` is the left-sidebar focus kind coordinated through ``ctx.focus_kind``) and
implement the hooks marked "subclass".
"""
from __future__ import annotations

import html

import numpy as np
from bokeh.models import BoxAnnotation, ColumnDataSource, Span, TextInput
from bokeh.plotting import figure

import abundance as ab
import store
from canvas import ACCENT, style_dark
from tableui import SortState, sort_headers
from tool import Tool
from units import to_range

import astropy.units as u

# scatter x-axis options: (key, axis label, per-line value fn) - shared by both wizards
X_AXES = {
    "evi": ("Excitation E_low [eV]", lambda ln: ln.e_low),
    "rew": ("Reduced EW  log(W/λ)", lambda ln: ab.reduced_ew(ln.ew, ln.lam0)),
    "lam0": ("Wavelength λ₀ [Å]", lambda ln: ln.lam0),
}

# a fitted profile (not trapz) is required - the FWHM window comes from it (see abundance.py)
PROFILES = ["gauss", "voigt", "asym_gauss"]
SYNTH_PAD_AA = 0.5           # extra Å each side of the auto window so the wings are synthesised
DEFAULT_JUMP_SPAN_AA = 1.0   # jump span for a no-window line when no line has been shown yet


def measure_candidate_lines(ref_spec, cand, *, profile, ew_min, snr, cont_pct):
    """Measure every candidate row of a master-linelist frame on the reference observation.

    The step-1 measuring loop shared by the abundance and stellar-parameter wizards:
    fit each line (:func:`abundance.measure_line`), keep no-core lines flagged but
    unmeasured, drop hard fit failures and clean lines below ``ew_min``, and measure the
    EW error (SNR + continuum-%) at selection — a line whose photon EW error exceeds its
    EW is flagged ``"ew_unc"`` and excluded. ``cand`` needs ``lam0``/``ion``/``evi``/
    ``loggf`` columns. Returns ``list[abundance.LineResult]``.
    """
    lines = []
    for row in cand.itertuples():
        line = {"lam0": float(row.lam0)}
        try:
            lprof, window, fit_window, flags = ab.measure_line(
                ref_spec, line, profile=profile)
        except Exception:
            continue          # hard fit failure (too few points) — drop it
        if "no_window" in flags or lprof is None:
            # no fittable core: list the line flagged, but measure no EW. Clicking it
            # just jumps to the line, keeping the previous view span (see _jump_line).
            lines.append(ab.LineResult(
                lam0=float(row.lam0), ion=str(row.ion), e_low=float(row.evi),
                loggf=float(row.loggf),
                mask=ab.LineMask(lam0=float(row.lam0), lmin=fit_window[0], lmax=fit_window[1],
                                 smin=window[0], smax=window[1], note=str(row.ion)),
                flags=tuple(flags), include=False))
            continue
        ew, fwhm, func = lprof.mAA, float(lprof.fwhm), lprof.func
        if not np.isfinite(ew):
            continue
        if not flags and ew < ew_min:
            continue          # an accepted line below the EW floor — drop as before
        # measure the EW error already at selection (SNR + continuum-%). A line whose
        # photon (SNR) EW error exceeds its own EW — fractional photon uncertainty > 1 —
        # is flagged and excluded (kept visible, like the other sanity rejections).
        # Keep the screened ``ew`` from measure_line; take only the errors here.
        ew_err = ew_err_phot = float("nan")
        try:
            prof = ab.measure_ew(ref_spec, window=fit_window, snr=snr,
                                 cont_pct=cont_pct, profile=profile, x0=float(row.lam0))
            ew_err, ew_err_phot = prof.width_err_mAA, prof.err_phot_mAA
        except Exception:
            pass
        flags = tuple(flags)
        if np.isfinite(ew_err_phot) and ew > 0 and ew_err_phot / ew > 1.0:
            flags = flags + ("ew_unc",)
        # sanity-check rejections (no fit window / off-centre fit / noisy) are kept but
        # excluded, so the user sees the flag and can override it in the table.
        lines.append(ab.LineResult(
            lam0=float(row.lam0), ion=str(row.ion), e_low=float(row.evi),
            loggf=float(row.loggf), ew=ew, ew_err=ew_err, fwhm=fwhm,
            mask=ab.LineMask(lam0=float(row.lam0), lmin=fit_window[0], lmax=fit_window[1],
                             smin=window[0], smax=window[1], note=str(row.ion)),
            fit_func=func, flags=flags, include=not flags))
    return lines


class AbundanceWizard(Tool):
    # --- subclass configuration --------------------------------------------------------
    prefix = ""            # DOM/bridge namespace: "abund" | "corr"
    kind = ""              # ctx.focus_kind value while a result owns the centre view
    btn_selector = ""      # toolbar button, e.g. "#abundance-btn"
    section = ""           # store section for the wizard's own form values
    empty_table_msg = "No result selected"
    y_axis_label = "A(X)"

    # sortable per-line columns shared by both wizards: {col: sort-value(ln)}. The header
    # markup lives in each subclass's table_html; here we own the state + the row order so
    # both tables sort identically (and stay in step with the Line Navigator). Wavelength /
    # excitation read low→high, EW and abundance magnitude are most useful largest-first.
    BASE_SORT_KEYS = {
        "lam0": lambda ln: ln.lam0,
        "ion": lambda ln: ln.ion,
        "e_low": lambda ln: ln.e_low,
        "loggf": lambda ln: ln.loggf,
        "ew": lambda ln: ln.ew,
        "a": lambda ln: ln.a,
    }
    BASE_SORT_ASC = {
        "lam0": True, "ion": True, "e_low": True, "loggf": False, "ew": False, "a": True,
    }
    # subclasses extend these for columns unique to their table (e.g. blend / vshift / vdop)
    extra_sort_keys = {}
    extra_sort_asc = {}

    def __init__(self, ctx, *, setup):
        Tool.__init__(self, ctx)
        self.sort_keys = {**self.BASE_SORT_KEYS, **self.extra_sort_keys}
        self.sort = SortState({**self.BASE_SORT_ASC, **self.extra_sort_asc}, "lam0")
        self.setup = setup                # SynthSetup: the wizard's own synthesis config
        self.own_mode_bridges((".sidebar-right", f"{self.prefix}-active"),
                              (self.btn_selector, "active"))
        self._view_bridge = ctx.bridges.cls("#main", f"{self.prefix}-mode")
        self._step_bridges = [ctx.bridges.cls("#main", f"{self.prefix}-step-{i}")
                              for i in (1, 2, 3)]
        self._results_bridge = ctx.bridges.html(f"{self.prefix}-results")
        self._table_bridge = ctx.bridges.html(f"{self.prefix}-table")

        self.step = 1
        self.show_toggled_only = False     # table header switch: hide excluded lines
        self._error = None                 # last worker error, surfaced on the doc thread
        self.results = []                  # list[ab.AbundanceResult]
        self.selected_uid = None
        self._last_span_aa = None          # span (Å) of the last shown line, for no-window jumps

        # step switch / run buttons / result clicks / line toggles all report through one
        # client->server channel, written as ``"<kind>:<arg>:<nonce>"``.
        self.actions = TextInput(value="", visible=False, name=f"{self.prefix}-actions")
        self.actions.on_change("value", self._on_action)

    # --- focus-kind coordination (which left-sidebar list owns the selection) -----------
    def is_focused(self):
        return getattr(self.ctx, "focus_kind", "spectrum") == self.kind

    def push_view(self, state):
        """Show/hide the wizard's centre view and hand the selection focus around."""
        from htmlbridge import pulse
        pulse(self._view_bridge, state)
        self.ctx.focus_kind = self.kind if state == "on" else "spectrum"
        self.ctx.emit("focus-changed", kind=self.ctx.focus_kind)

    def set_step(self, n):
        from htmlbridge import toggle
        self.step = max(1, min(3, int(n)))
        for i, bridge in enumerate(self._step_bridges, start=1):
            toggle(bridge, i == self.step)

    # --- dark scatter figure (centre view, step 3) ---------------------------------------
    def _build_fig(self, source):
        self.src = source
        fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                     tools="pan,wheel_zoom,box_zoom,reset", active_scroll="wheel_zoom",
                     min_border_left=62, min_border_bottom=42)
        style_dark(fig)   # match the other Bokeh panels' dark theme
        fig.scatter("x", "a", source=self.src, size=8, fill_color=ACCENT,
                    line_color="#fff", line_width=1)
        self.mean_line = Span(location=0, dimension="width", line_color="#e05c5c",
                              line_width=1.2, line_dash="dashed", visible=False)
        self.mean_band = BoxAnnotation(fill_color="#e05c5c", fill_alpha=0.08, visible=False)
        fig.add_layout(self.mean_line)
        fig.add_layout(self.mean_band)
        fig.yaxis.axis_label = self.y_axis_label
        fig.xaxis.axis_label = X_AXES[str(self.form3.value("xaxis"))][0]
        fig.name = f"{self.prefix}_fig"
        self.fig = fig
        return fig

    # --- Tool protocol -------------------------------------------------------------------
    @property
    def models(self):
        return [*self.form1.models, *self.form3.models, self.actions, *self.setup.models]

    def activate(self):
        # first open: seed the wizard's own synthesis setup from the Synthesize tab, then
        # stay independent. Detect first open with a marker (not section-emptiness) because
        # the SynthSetup constructor pre-writes a validated atmos_grid_id.
        if not store.get_section(self.setup.section).get("_seeded"):
            self.setup.inherit_from()
            store.update_pref(self.setup.section, "_seeded", True)
        self.active = True
        self._on_activate_extra()
        self.setup.render()               # step-2 forms + grid + atoms note
        self.push_mode("on")
        self.set_step(self.step)
        # if a result is already selected, re-enter its centre view
        if self.selected() is not None and self.is_focused():
            self.push_view("on")

    def _on_activate_extra(self):
        """Subclass hook: extra surfaces to push on activate (e.g. the reference selector)."""

    def deactivate(self):
        was = self.active
        self.active = False
        if was:
            self.push_view("off")
            self._clear_line_decor()
            self.push_mode("off")

    def _clear_line_decor(self):
        self.canvas.clear_abund_overlay()
        self.canvas.clear_abund_hatch()
        self.canvas.clear_abund_window()

    # --- the shared slice of the action channel ------------------------------------------
    def _on_action(self, attr, old, new):
        if not new:
            return
        parts = new.split(":")             # "<kind>:<arg>:<nonce>"
        kind = parts[0]
        arg = parts[1] if len(parts) > 1 else ""
        if kind == "step":
            self.set_step(int(arg))
        elif kind == "select":
            self.select_result(arg)
        elif kind == "remove":
            self.remove_result(arg)
        elif kind == "toggle-line":
            self._toggle_line(arg)
        elif kind == "mask-toggle":
            self.show_toggled_only = not self.show_toggled_only
            self._render_selected()
        elif kind == "jump-line":
            self._jump_line(arg)
        elif kind == "sort":
            if self.sort.toggle(arg):
                res = self.selected()
                if res is not None:
                    self.push_table(self.table_html(res))
        else:
            self._extra_action(kind, arg)

    def _extra_action(self, kind, arg):
        """Subclass hook: the wizard-specific run actions."""

    # --- results list / selection ---------------------------------------------------------
    def _find(self, uid):
        return next((r for r in self.results if r.uid == uid), None)

    def selected(self):
        return self._find(self.selected_uid)

    def has_results(self):
        return bool(self.results)

    def select_result(self, uid):
        if self._find(uid) is None:
            return
        self.selected_uid = uid
        self.ctx.tools.activate(self.name)   # a result click opens the panel (no toolbar needed)
        self._clear_line_decor()
        self.push_results()
        self._render_selected()
        self.push_view("on")

    def blur(self):
        """A spectrum object took focus: leave the centre view but keep the selection."""
        self.push_view("off")
        self._clear_line_decor()
        self.push_results()                  # drop the result highlight

    def remove_result(self, uid):
        was_selected = self.selected_uid == uid
        self.results = [r for r in self.results if r.uid != uid]
        if was_selected:
            self.selected_uid = self.results[-1].uid if self.results else None
        self.push_results()
        if self.selected_uid and was_selected:
            self.select_result(self.selected_uid)
        elif not self.results:
            self.push_view("off")
            self._clear_line_decor()

    def _toggle_line(self, token):
        res = self.selected()
        if res is None:
            return
        try:
            idx = int(token)
        except ValueError:
            return
        if 0 <= idx < len(res.lines):
            res.lines[idx].include = not res.lines[idx].include
            if res.stage == self.final_stage:
                res.summary = self._summarise(res.lines)
                res.name = self._result_label(res)
                self.push_results()
            self._render_selected()

    # subclass: the stage name after the last run, the aggregation and the sidebar label
    final_stage = "uncert"

    def _summarise(self, lines):
        raise NotImplementedError

    def _result_label(self, res):
        raise NotImplementedError

    # --- streaming + rendering -------------------------------------------------------------
    def _stream_row(self, res):
        """Doc-thread callback: repaint the table as each line's synthesis lands."""
        if self.selected() is res and self.is_focused():
            self.push_table(self.table_html(res))

    def push_results(self):
        self._results_bridge.value = self.results_html()

    def push_table(self, html_text):
        self._table_bridge.value = html_text

    # --- shared sortable-table scaffolding ------------------------------------------------
    def _line_order(self, res):
        """Indices into ``res.lines`` in the active sort order (NaN/None keys sort last).

        Rows carry their *original* index as ``data-line`` so the client's toggle/jump
        actions keep addressing ``res.lines[i]`` no matter how the display is sorted.
        """
        return self.sort.order(res.lines, self.sort_keys[self.sort.col])

    def _table_head(self, mask_class, columns):
        """The clickable ``<thead>`` row: a leading mask checkbox, then ``columns``.

        ``columns`` is a sequence of ``(label, sort-col | None)`` pairs; a ``None`` column
        (e.g. the flags column) renders as a plain, non-sortable header.
        """
        checked = " checked" if self.show_toggled_only else ""
        chk = (f'<input type="checkbox" class="{mask_class}"{checked} '
               'title="Show only included lines">')
        return "<tr>" + sort_headers(self.sort, [(chk, None), *columns]) + "</tr>"

    def _render_selected(self):
        res = self.selected()
        if res is None:
            self._render_scatter()
            self.push_table(f'<div class="linelist-empty">{self.empty_table_msg}</div>')
            return
        self._render_scatter()
        self.push_table(self.table_html(res))

    def _crop_scatter(self, xs, ys):
        """Fix the scatter ranges to the shown values with a 5% margin each side."""
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        xs, ys = xs[np.isfinite(xs)], ys[np.isfinite(ys)]
        if xs.size == 0 or ys.size == 0:
            return
        for rng, lo, hi in ((self.fig.x_range, xs.min(), xs.max()),
                            (self.fig.y_range, ys.min(), ys.max())):
            pad = 0.05 * (hi - lo) or 0.05 * (abs(hi) or 1.0)
            rng.start, rng.end = float(lo - pad), float(hi + pad)

    # --- jump-to-line -----------------------------------------------------------------------
    def _jump_line(self, token):
        res = self.selected()
        if res is None:
            return
        try:
            ln = res.lines[int(token)]
        except (ValueError, IndexError):
            return
        lo, hi = ln.window
        if not np.isfinite(lo):
            # no window was determined for this line - just recentre on it, keeping the span
            # from the previously shown line, and drop the (absent) fit overlays.
            span = self._last_span_aa or DEFAULT_JUMP_SPAN_AA
            self._set_range(ln.lam0 - 0.5 * span, ln.lam0 + 0.5 * span)
            self._clear_line_decor()
            return
        self._last_span_aa = abs(hi - lo)
        self._set_range(lo, hi)
        self._jump_decor(res, ln)

    def _jump_decor(self, res, ln):
        """Subclass hook: the overlays/hatch/window bands drawn for a jumped-to line."""

    def _set_range(self, lo_aa, hi_aa):
        """Set the main + nav x-range to the Å interval (display units, clamped to bounds)."""
        a, b = to_range(lo_aa, hi_aa, u.AA, self.canvas.xunit)
        lo_d, hi_d = min(a, b), max(a, b)
        for rng in (self.canvas.main.x_range, self.canvas.nav.x_range):
            s, e = lo_d, hi_d
            if rng.bounds:
                blo, bhi = min(rng.bounds), max(rng.bounds)
                s, e = max(s, blo), min(e, bhi)
            if s < e:
                rng.start, rng.end = s, e

    # --- shared formatting --------------------------------------------------------------------
    @staticmethod
    def _fmt(value, err=None, plus=False, digits=3):
        """Value ± err, or an em-dash when the stage that fills it hasn't run yet."""
        if value is None or not np.isfinite(value):
            return "—"
        sign = "+" if plus else ""
        out = f"{value:{sign}.{digits}f}"
        if err is not None and np.isfinite(err):
            out += f" ± {err:.{digits}f}"
        return out

    def _results_rows(self, badge, item_cls, data_key, remove_key):
        focused = self.is_focused()
        rows = []
        for r in self.results:
            sel = " selected" if (r.uid == self.selected_uid and focused) else ""
            rows.append(
                f'<div class="object-item result-item{item_cls}{sel}" {data_key}="{r.uid}">'
                f'<span class="result-badge">{badge}</span>'
                f'<span class="object-name" title="{html.escape(r.name)}">'
                f'{html.escape(r.name)}</span>'
                f'<button class="result-btn" {remove_key}="{r.uid}" '
                f'title="Remove"><i class="ti ti-x"></i></button></div>')
        return rows
