"""The shared machinery of the stepped abundance wizards.

:class:`abundancetool.AbundanceTool` (element abundances vs an observation) and
:class:`abundcorrtool.AbundCorrTool` (NLTE-LTE corrections vs a synthetic 1D LTE
reference) grew as siblings with near-verbatim copies of the step switcher, the
results-list / selection management, the action-channel dispatch, the per-line table
scaffolding, the dark scatter figure and the jump-to-line logic. That common machinery
lives here; the subclasses keep what is genuinely different - their forms, how many steps
they have, their run stages (what gets synthesised and inverted) and their table columns.

Subclasses set the class attributes below (``prefix`` namespaces the DOM/bridge ids,
``kind`` is the left-sidebar focus kind coordinated through ``ctx.focus_kind``) and
implement the hooks marked "subclass".
"""
from __future__ import annotations

import html
from pathlib import Path
from dataclasses import replace

import numpy as np
from bokeh.models import (BoxAnnotation, ColumnDataSource, CustomJS, HoverTool, Range1d,
                          Span, TextInput, Whisker, PointDrawTool)
from bokeh.plotting import figure

from tumult import abundance as ab
from tumult.abundance import SYNTH_PAD_AA  # noqa: F401 (the tools import it from here)
import store
from canvas import ACCENT, style_dark
from tableui import SortState, sort_headers
from tool import Tool
from units import to_range, to_value

import astropy.units as u

# scatter x axes: key -> (axis label, per-line value fn)
X_AXES = {
    "evi": ("Excitation E_low [eV]", lambda ln: ln.e_low),
    "rew": ("Reduced EW  log(W/λ)", lambda ln: ab.reduced_ew(ln.ew, ln.lam0)),
}

# the two diagnostics every abundance fit is read against, shown as a pair of panels: the
# abundance against excitation energy (a temperature diagnostic) and against reduced EW (a
# microturbulence one).
PANEL_AXES = ("evi", "rew")

# a fitted profile (not trapz) is required - the FWHM window comes from it (see abundance.py)
PROFILES = ["gauss", "voigt", "asym_gauss"]
DEFAULT_JUMP_SPAN_AA = 1.0   # jump span for a no-window line when no line has been shown yet
FIT_WINDOW_MARGIN = 0.05     # fraction of a fitted line window shown beyond either edge

# The "Local continuum" header (abund/ionbal only) is one of the table's non-sortable
# columns, so the shared line-list rule would otherwise left-align it while its checkboxes
# below stay right-aligned (see .local-continuum-head in style.css). A class on the label
# itself, not a column-position selector, so it survives a column being added or reordered.
LOCAL_CONTINUUM_HEAD = '<span class="local-continuum-head">Local continuum</span>'

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

# The shared client->server wiring of every wizard and of the ionisation-balance tool,
# emitted once through the ``analysis_script`` template variable. Each tool's markup scopes
# itself with ``data-tool="<prefix>"`` (the step tabs, the results rows, the lines table);
# the binder reads the prefix from the closest scope and writes ``"<kind>:<arg>:<nonce>"``
# to the ``<prefix>-actions`` channel, which :meth:`AbundanceWizard._on_action` splits on
# colons (uids and line indices contain none). A row click jumps to its line; a click on the
# row's include checkbox only toggles the line. The tools' own run/refit buttons and the
# atom sidebar (which removes its nested corrections) stay in the page script.
ANALYSIS_SCRIPT = r"""
(function () {
  var seq = 0;   // a repeated identical action within one millisecond is still a new value
  function analysisAction(prefix, kind, arg) {
    var m = Bokeh.documents[0].get_model_by_name(prefix + '-actions');
    if (m) m.setv({ value: kind + ':' + (arg || '') + ':' + Date.now() + '.' + (++seq) });
  }
  document.addEventListener('click', function (ev) {
    var scope = ev.target.closest('[data-tool]');
    if (!scope) return;
    var prefix = scope.dataset.tool;
    var step = ev.target.closest('.step-tab[data-step]');
    if (step) { analysisAction(prefix, 'step', step.dataset.step); return; }
    var rem = ev.target.closest('[data-result-remove]');
    if (rem) { analysisAction(prefix, 'remove', rem.dataset.resultRemove); return; }
    var item = ev.target.closest('[data-result]');
    if (item) { analysisAction(prefix, 'select', item.dataset.result); return; }
    var sort = ev.target.closest('.line-sort[data-sortcol]');
    if (sort) { analysisAction(prefix, 'sort', sort.dataset.sortcol); return; }
    var remaining = ev.target.closest('[data-complete-rungs]');
    if (remaining) { analysisAction(prefix, 'complete-rungs'); return; }
    var refit = ev.target.closest('[data-refit-continuum]');
    if (refit) { analysisAction(prefix, 'refit-continuum'); return; }
    var row = ev.target.closest('tr[data-line]');
    if (!row || ev.target.closest('.line-include-chk, .line-continuum-chk')) return;
    scope.querySelectorAll('tr.row-selected')
         .forEach(function (r) { r.classList.remove('row-selected'); });
    row.classList.add('row-selected');
    analysisAction(prefix, 'jump-line', row.dataset.line);
  });
  document.addEventListener('change', function (ev) {
    var scope = ev.target.closest('[data-tool]');
    if (!scope) return;
    var prefix = scope.dataset.tool;
    if (ev.target.closest('.analysis-mask-chk')) { analysisAction(prefix, 'mask-toggle'); return; }
    var cont = ev.target.closest('.line-continuum-chk');
    if (cont) { analysisAction(prefix, 'continuum', cont.dataset.line); return; }
    var chk = ev.target.closest('.line-include-chk');
    if (chk) analysisAction(prefix, 'toggle-line', chk.dataset.line);
  });
})();
"""


from jobs import ResultJobs


def preparation_note(res):
    """The processing a result's lines were measured under (``MeasuredLines.preparation``,
    written by the measuring tool from ``canvas.preparation``), as a hint; empty when the
    result predates that record."""
    prep = getattr(res, "preparation", None) or {}
    if not prep:
        return ""
    parts = [str(prep.get("normalisation", "")).capitalize(),
             "filtered" if prep.get("filtered") else "unfiltered"]
    if prep.get("line_list"):
        parts.append(f"line list {Path(str(prep['line_list'])).name}")
    return f'<div class="hint">Measured with: {html.escape(" · ".join(p for p in parts if p))}.</div>'


class AbundanceWizard(Tool):
    # --- subclass configuration --------------------------------------------------------
    prefix = ""            # DOM/bridge namespace: "abund" | "corr"
    steps = 3              # how many step panels the wizard's switcher drives
    kind = ""              # ctx.focus_kind value while a result owns the centre view
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
        self.own_mode_bridges((".sidebar-right", f"{self.prefix}-active"))
        self._view_bridge = ctx.bridges.cls("#main", f"{self.prefix}-mode")
        self._step_bridges = [ctx.bridges.cls("#main", f"{self.prefix}-step-{i}")
                              for i in range(1, self.steps + 1)]
        self._results_bridge = ctx.bridges.html(f"{self.prefix}-results")
        self._table_bridge = ctx.bridges.html(f"{self.prefix}-table")
        # another list took (or gave back) the single left-sidebar selection: restyle our rows
        ctx.subscribe("focus-changed", lambda **_: self.push_results())

        self.step = 1
        self.show_toggled_only = False     # table header switch: hide excluded lines
        self.results = []                  # list[ab.AbundanceResult]
        self.jobs = ResultJobs(ctx, self._find, self._job_progress)
        if self.task_kind:
            # a record of this kind another document submitted lands here on a reconnect
            ctx.task_importers[self.task_kind] = self.import_task_result
            ctx.task_progress[self.task_kind] = self.import_task_progress
        self.selected_uid = None
        self._last_span_aa = None          # span (Å) of the last shown line, for no-window jumps
        self._jumped = None                # index of the line whose window is accented

        # step switch / run buttons / result clicks / line toggles all report through one
        # client->server channel, written as ``"<kind>:<arg>:<nonce>"``.
        self.actions = TextInput(value="", visible=False, name=f"{self.prefix}-actions")
        self.actions.on_change("value", self._on_action)
        self._continuum_sync = False
        self.continuum_source = ColumnDataSource(dict(x=[], y=[]))
        # These are editable markers, not a spectrum: canvas.py keeps them off the main
        # figure's y auto-fit (see its note on ``main.y_range.renderers``), so dragging one
        # cannot widen the range to chase the very point being moved.
        self.continuum_markers = self.canvas.main.scatter(
            "x", "y", source=self.continuum_source, size=10, color=ACCENT, visible=False)
        self.continuum_curve = self.canvas.main.line(
            "x", "y", source=self.continuum_source, color=ACCENT,
            line_dash="dashed", visible=False)
        self.continuum_drag = PointDrawTool(renderers=[self.continuum_markers], add=False)
        self.canvas.main.add_tools(self.continuum_drag)
        self.continuum_source.on_change("data", self._on_continuum_drag)

    # --- window highlights ---------------------------------------------------------------
    # A result's windows belong to the result, not to a jump: a wizard that paints them shows
    # them for as long as it owns the centre view, with the jumped-to line accented. Only the
    # corrections wizard lets that accent be dragged (its fit windows are its own to edit);
    # the others' windows are edited in the LineMasks view. The call sites are the
    # subclass's, so a wizard whose lines carry no user-drawn windows paints nothing.
    def _result_masks(self, res):
        """The masks to highlight for ``res``: its lines' windows, carrying the *line's*
        include state, so a row switched off in the table is drawn switched off."""
        return [replace(ln.mask, include=ln.include) for ln in res.lines]

    def _result_spec(self, res):
        """Subclass hook: the spectrum ``res`` was measured on, whose points snap the window
        bounds, or None when it is no longer loaded."""
        return None

    def _render_windows(self, active=None):
        res = self.selected()
        if not (self.active and self.is_focused()) or res is None:
            return self.canvas.clear_windows()
        spec = self._result_spec(res)
        xpts = np.asarray(spec.values[0], dtype=float) if spec is not None else None
        self.canvas.highlight_windows(self._result_masks(res), active=active, xpts=xpts)

    # --- synthetic overlays ----------------------------------------------------------------
    # Like the windows, the syntheses belong to the result: a wizard that keeps one per line
    # draws them all for as long as it owns the centre view, so the spectrum pane shows what
    # the fit is standing on and fills in as the run lands. Jumping a row then only moves the
    # accent and the EW hatch, it does not decide what is drawn.
    paint_overlays = False   # subclass: True where every line carries a synthesis worth drawing

    def _render_decor(self, active=None):
        """Everything the selected result puts on the spectrum: its windows (``active``
        accented) and the syntheses standing on them."""
        self._render_windows(active)
        self._render_overlays()
        self._render_continuum()

    def _render_overlays(self, res=None):
        if not self.paint_overlays:
            return
        res = res if res is not None else self.selected()
        if not (self.active and self.is_focused()) or res is None:
            return self.canvas.clear_abund_overlay()
        self.canvas.set_abund_overlays([(ab.displayed_model(ln), ln.window) for ln in res.lines
                                        if ln.include and ln.fit_spec is not None])

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
        self.step = max(1, min(self.steps, int(n)))
        for i, bridge in enumerate(self._step_bridges, start=1):
            toggle(bridge, i == self.step)

    # --- dark scatter figures (the centre view's diagnostic pane) ------------------------
    # Every wizard shows the pair of diagnostics (:data:`PANEL_AXES`) side by side
    # (:meth:`_build_figs`), drawn from one ColumnDataSource carrying an ``x_<axis>``
    # column per panel, so a single data push repaints both.
    def _new_fig(self, xkey):
        fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                     tools="pan,wheel_zoom,box_zoom,tap,reset", active_scroll="wheel_zoom",
                     min_border_left=62, min_border_bottom=42,
                     x_range=Range1d(0.0, 1.0), y_range=Range1d(0.0, 1.0))
        style_dark(fig)   # match the other Bokeh panels' dark theme
        span = Span(location=0, dimension="width", line_color="#e05c5c",
                    line_width=1.2, line_dash="dashed", visible=False)
        band = BoxAnnotation(fill_color="#e05c5c", fill_alpha=0.08, visible=False)
        fig.add_layout(span)
        fig.add_layout(band)
        fig.yaxis.axis_label = self.y_axis_label
        self.figs.append(fig)
        self.mean_lines.append(span)
        self.mean_bands.append(band)
        return fig

    def _build_figs(self, source):
        """The two diagnostic panels (:data:`PANEL_AXES`), side by side."""
        self.src = source
        source.name = f"{self.prefix}_points"
        self.figs, self.mean_lines, self.mean_bands = [], [], []
        for i, key in enumerate(PANEL_AXES):
            fig = self._new_fig(f"x_{key}")
            fig.xaxis.axis_label = X_AXES[key][0]
            fig.name = f"{self.prefix}_fig" + ("" if i == 0 else str(i + 1))
        # the ring the page script moves onto the point of a hovered table row, on every
        # panel at once (it carries both panels' x, like the points themselves)
        self.hover_src = ColumnDataSource({f"x_{k}": [] for k in PANEL_AXES} | dict(a=[]))
        self.hover_src.name = f"{self.prefix}_hover"
        for key, fig in zip(PANEL_AXES, self.figs):
            fig.scatter(f"x_{key}", "a", source=self.hover_src, size=16, fill_color=None,
                        line_color="#fff", line_width=2)
        self._add_series(source, ACCENT)
        self.fig, self.fig2 = self.figs
        return self.figs

    def _add_series(self, source, color, xkeys=None):
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

    def _add_whiskers(self, source, color):
        """Error bars on both panels, from ``source``'s ``lo``/``hi`` columns."""
        for key, fig in zip(PANEL_AXES, self.figs):
            fig.add_layout(Whisker(base=f"x_{key}", upper="hi", lower="lo", source=source,
                                   line_color=color, line_width=1.2))

    def _pack_points(self, indexed, value_of=lambda ln: ln.a, sigma_of=lambda ln: ln.sigma_a):
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

    @staticmethod
    def empty_points():
        return ({f"x_{k}": [] for k in PANEL_AXES}
                | dict(a=[], lo=[], hi=[], lam=[], ion=[], idx=[]))

    @staticmethod
    def value_label(res):
        """What a result's ``ln.a`` actually is, for the table header, the axis and the mean.

        A run whose abundances were given relative to the model composition inverts each line
        to [X/Fe] (``A(X) = A_solar(X) + [Fe/H] + [X/Fe]``), not to an absolute abundance;
        for iron itself that value is the offset from the model's own [Fe/H], so it is named
        Δ[Fe/H] rather than the meaningless [Fe/Fe]. An absolute run (the balance always,
        which converts every inverted line back to A(X)) is labelled A(X).
        """
        el = getattr(res, "symbol", "") or getattr(res, "element", "") or "X"
        if not getattr(res, "relative", False):
            return f"A({el})"
        return "Δ[Fe/H]" if el == "Fe" else f"[{el}/Fe]"

    def _crop_panels(self, series):
        """Fit every panel to its x data and the full vertical whisker span.

        Bokeh's automatic ranges do not account for ``Whisker`` annotations.  Derive the
        shared y extent here from the lower and upper error-bar columns so callers cannot
        accidentally crop to the point centres and clip a bar.
        """
        ys = [y for data in series for column in ("lo", "hi") for y in data[column]]
        for key, fig in zip(PANEL_AXES, self.figs):
            xs = [x for d in series for x in d[f"x_{key}"]]
            self._crop_scatter(xs, ys, fig)

    def _set_mean(self, value, err=0.0):
        """Show the mean (and its ± band) on every panel, or hide it when there is none."""
        ok = value is not None and np.isfinite(value)
        err = err if np.isfinite(err) else 0.0
        for span, band in zip(self.mean_lines, self.mean_bands):
            span.visible = band.visible = bool(ok)
            if ok:
                span.location = float(value)
                band.top, band.bottom = float(value) + err, float(value) - err

    # --- Tool protocol -------------------------------------------------------------------
    @property
    def models(self):
        return [*self.form1.models, self.actions, *self.setup.models]

    def activate(self):
        # first open: seed the wizard's own synthesis setup from the synthesis modal, then
        # stay independent. Detect first open with a marker (not section-emptiness) because
        # the SynthSetup constructor pre-writes a validated atmos_grid_id.
        if not store.get_section(self.setup.section).get("_seeded"):
            self.setup.inherit_from()
            store.update_pref(self.setup.section, "_seeded", True)
        self.active = True
        self.setup.render()               # the synthesis forms + grid + atoms note
        self.push_mode("on")
        self.set_step(self.step)
        # if a result is already selected, re-enter its centre view
        if self.selected() is not None and self.is_focused():
            self.push_view("on")

    def deactivate(self):
        was = self.active
        self.active = False
        if was:
            self.push_view("off")
            self._clear_line_decor()
            self.push_mode("off")

    def _clear_line_decor(self):
        self._hide_continuum()
        self.canvas.clear_abund_overlay()
        self.canvas.clear_abund_hatch()
        self.canvas.clear_windows()

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
        elif kind == "continuum":
            self._toggle_continuum(arg)
        elif kind == "complete-rungs":
            self._complete_rungs()
        elif kind == "refit-continuum":
            if self.prefix == "abund":
                self.run_synth(request=self.selected().request)
            elif self.prefix == "ionbal":
                self.open_refit()
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

    # --- the shared EW refresh --------------------------------------------------------------
    def report_refresh(self, res, outcome):
        """Repaint ``res`` after :func:`tumult.abundance.refresh_ews` changed any of its
        lines, and name the lines whose measurement failed (they keep their old EW)."""
        if outcome.changed:
            self._render_selected()
        if outcome.failed:
            lams = ", ".join(f"{res.lines[i].lam0:.2f}" for i in sorted(outcome.failed))
            self.notify(f"{res.name}: EW refresh failed for {len(outcome.failed)} "
                        f"line{'s' if len(outcome.failed) != 1 else ''} ({lams} Å); the "
                        "previous values were kept.")

    # --- results list / selection ---------------------------------------------------------
    def _find(self, uid):
        return next((r for r in self.results if r.uid == uid), None)

    # --- durable calculations: this wizard's task kind -----------------------------------
    # The wizard's calculation runs as a workflow record of ``task_kind``
    # (tumult.synthesis.tasks), so it outlives the tab and the server. Its snapshots
    # and result are result objects carrying the uid of the one they continue; a uid this
    # list does not hold (a record another document submitted) is adopted as a new
    # result. Subclasses that submit through ``self.jobs.submit_operation`` name the operation.
    task_kind = None

    def _adopt(self, incoming):
        """The listed result ``incoming`` continues, updated in place from it, or
        ``incoming`` itself once listed and selected."""
        res = self._find(incoming.uid)
        if res is None:
            self.results.append(incoming)
            self.selected_uid = incoming.uid
            self.push_results()
            return incoming
        res.__dict__.update(incoming.__dict__)
        return res

    def import_task_progress(self, snapshot, handle):
        self._job_progress(self._adopt(snapshot))

    def import_task_result(self, result, handle):
        res = self._adopt(result)
        res.name = self._result_label(res)
        self._job_progress(res)

    def selected(self):
        return self._find(self.selected_uid)

    # --- session persistence: the envelope every result list shares ----------------------
    # The result objects themselves are built by the subclass's factory and relinked by
    # its ``relink``; only the list/selected-uid envelope and its validation live here.
    result_from_dict = staticmethod(ab.AbundanceResult.from_dict)

    def snapshot(self):
        """``{"results": [...], "selected": uid}`` for the session slot."""
        return {"results": [r.to_dict() for r in self.results], "selected": self.selected_uid}

    def deserialise(self, results):
        """The result objects of a saved envelope (no live state touched); a malformed
        record raises."""
        return [self.result_from_dict(d) for d in results]

    def relink(self, results, spectrum_uids, maskset_uids, atom_paths):
        """Clear the references of ``results`` that name objects a load does not restore;
        returns one line per cleared reference. The base clears the reference observation
        (``ref_uid``); subclasses add their own links."""
        notes = []
        for res in results:
            if res.ref_uid is not None and res.ref_uid not in spectrum_uids:
                notes.append(f"{self.name} result {res.name or res.uid}: reference spectrum "
                             f"{res.ref_name or res.ref_uid} is not loaded")
                res.ref_uid = None
        return notes

    def restore(self, results, selected):
        """Replace the result list (and repaint the sidebar); ``selected`` must be one of
        them or ``None``. Pending result jobs of the old list are forgotten."""
        self.jobs.current.clear()
        self.results = list(results)
        uids = {r.uid for r in self.results}
        self.selected_uid = selected if selected in uids else None
        self._jumped = None
        self.push_results()

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
        self.jobs.invalidate(uid)
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
            self.jobs.invalidate(res.uid)
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
    def _job_progress(self, res):
        self.push_results()
        self._stream_row(res)

    def _stream_row(self, res):
        """Doc-thread callback: repaint everything an incoming line changed — the table, the
        diagnostic panels and the synthetic overlays — so a fit is watched as it runs."""
        if self.selected() is res and self.is_focused():
            self.push_table(self.table_html(res))
            self._render_scatter()
            self._render_overlays(res)

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

    def _table_head(self, columns):
        """The clickable ``<thead>`` row: a leading mask checkbox, then ``columns``.

        ``columns`` is a sequence of ``(label, sort-col | None)`` pairs; a ``None`` column
        (e.g. the flags column) renders as a plain, non-sortable header.
        """
        checked = " checked" if self.show_toggled_only else ""
        chk = (f'<input type="checkbox" class="analysis-mask-chk"{checked} '
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

    def _crop_scatter(self, xs, ys, fig=None):
        """Fix a panel's ranges to the shown values with a 5% margin each side."""
        fig = fig if fig is not None else self.fig
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
        lo, hi = min(lo, hi), max(lo, hi)
        self._last_span_aa = hi - lo
        pad = FIT_WINDOW_MARGIN * (hi - lo)
        self._set_range(lo - pad, hi + pad)
        self._jump_decor(res, ln)

    def _jump_decor(self, res, ln):
        """Subclass hook: the overlays/hatch/window highlights drawn for a jumped-to line."""

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

    def _results_rows(self, badge, results=None, nested=False):
        """Sidebar rows for ``results`` (default: all); ``nested`` indents them under the
        object they were computed from. The row carries ``data-tool`` (this wizard's
        prefix) and ``data-result`` (the uid), its button ``data-result-remove``: the
        attributes the page's analysis binder routes select/remove clicks by."""
        focused = self.is_focused()
        rows = []
        nest = " result-nested" if nested else ""
        for r in (self.results if results is None else results):
            sel = " selected" if (r.uid == self.selected_uid and focused) else ""
            rows.append(
                f'<div class="object-item result-item{nest}{sel}" '
                f'data-tool="{self.prefix}" data-result="{r.uid}">'
                f'<span class="result-badge">{badge}</span>'
                f'<span class="object-name" title="{html.escape(r.name)}">'
                f'{html.escape(r.name)}</span>'
                f'<button class="result-btn" data-result-remove="{r.uid}" '
                f'title="Remove from workspace"><i class="ti ti-x"></i></button></div>')
        return rows

    # --- measurement-owned local continua ------------------------------------------
    def _hide_continuum(self):
        self.continuum_markers.visible = self.continuum_curve.visible = False
        self.canvas.release_edit_tool(self.continuum_drag)

    def _render_continuum(self):
        res = self.selected()
        if (self.prefix not in ("abund", "ionbal") or not self.active
                or not self.is_focused() or res is None or self._jumped is None
                or self._jumped >= len(res.lines)):
            return self._hide_continuum()
        ln = res.lines[self._jumped]
        if not ln.renorm or not ln.continuum:
            return self._hide_continuum()
        xx, yy = np.asarray(ln.continuum).T
        self._continuum_sync = True
        try:
            self.continuum_source.data = dict(x=to_value(xx, u.AA, self.canvas.xunit), y=yy)
        finally:
            self._continuum_sync = False
        self.continuum_markers.visible = self.continuum_curve.visible = True
        self.canvas.claim_edit_tool(self.continuum_drag)
        self.canvas.set_abund_hatch(ln.fit_func, ln.fit_window, nodes=ln.continuum)

    def _continuum_changed(self, res, ln, **kwargs):
        spec = self._result_spec(res)
        if spec is None:
            self.notify("The target spectrum is no longer loaded.")
            return
        try:
            ab.set_local_continuum(res, ln, spec, **kwargs)
        except (ValueError, RuntimeError) as exc:
            self.notify(str(exc))
            self._render_continuum()
            return
        self.jobs.invalidate(res.uid)
        res.name = self._result_label(res)
        self.push_results()
        self._render_selected()
        self._jump_decor(res, ln)

    def _toggle_continuum(self, token):
        res = self.selected()
        if res is None or self.prefix not in ("abund", "ionbal"):
            return
        idx = int(token)
        if not 0 <= idx < len(res.lines):
            return
        self._jumped = idx
        ln = res.lines[idx]
        self._continuum_changed(res, ln, enabled=not ln.renorm)
        self._jump_line(token)

    def _on_continuum_drag(self, attr, old, new):
        if self._continuum_sync or not self.continuum_markers.visible:
            return
        res = self.selected()
        if res is None or self._jumped is None:
            return
        ln = res.lines[self._jumped]
        if len(new.get("y", [])) != 2:
            self._render_continuum()
            return
        # The markers move vertically; their wavelengths remain at the synth edges.
        nodes = tuple(zip(ln.window, (float(y) for y in new["y"])))
        self._continuum_changed(res, ln, nodes=nodes)

    @staticmethod
    def continuum_cell(ln, index):
        checked = " checked" if ln.renorm else ""
        return (f'<td><input type="checkbox" class="line-continuum-chk" '
                f'data-line="{index}"{checked} title="Adjust this measurement’s local continuum">'
                '</td>')

    def continuum_actions(self, res):
        if any(ln.missing_rungs for ln in res.lines):
            return (f'<div data-tool="{self.prefix}"><span class="hint">'
                    'The edited EW needs additional abundance rungs. </span>'
                    '<button class="primary-btn" data-complete-rungs>'
                    'Compute remaining abundance rungs</button></div>')
        if any("refit_required" in ln.flags for ln in res.lines):
            return (f'<div data-tool="{self.prefix}"><button class="primary-btn" '
                    'data-refit-continuum>Refit with edited continuum</button></div>')
        return ""

    def _complete_rungs(self):
        res = self.selected()
        if res is None or not any(ln.missing_rungs for ln in res.lines):
            return
        from copy import deepcopy
        target = deepcopy(self._result_spec(res))
        def work(working, token, progress):
            return ab.complete_result_rungs(working, target=target,
                                            cancel_token=token, progress=progress)
        def done(result):
            result.name = self._result_label(result)
            self.push_results()
            self._render_selected()
            self.notify("Cached curves of growth updated." if not any(
                ln.missing_rungs for ln in result.lines) else
                "More abundance rungs are still needed; inspect the flagged lines.")
        self.jobs.submit("Complete abundance rungs", work, done, result=res, detached=True)
