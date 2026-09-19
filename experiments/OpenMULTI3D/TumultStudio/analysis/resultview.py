"""The result view every analysis tool presents: navigation, the per-line table and the
spectrum decoration that belongs to the result in view.

:class:`abundancetool.AbundanceTool` (element abundances against an observation),
:class:`abundcorrtool.AbundCorrTool` (NLTE−LTE corrections against a synthetic 1D LTE
reference) and :class:`ionbalance.IonBalanceTool` (stellar parameters from an ionisation
balance) each produce a different scientific object, but a user reads all three the same
way: pick a result in the left sidebar, step between its panels, sort and curate its lines,
click a row to jump the spectrum to that line, and see whether the solution still stands on
its inputs. :class:`ResultView` owns that: the step switcher, the single client→server
action channel (:data:`ANALYSIS_SCRIPT`), the selection and identity operations over its
:class:`~.results.ResultCollection`, the sortable table scaffolding with its checkbox
columns, the jump-to-line gesture and the out-of-date banner.

It composes rather than contains the rest: :class:`~.diagnostics.DiagnosticPanels` for the
two scatter panels, :class:`~.localcontinuum.LocalContinuum` for the per-line continuum
(the views whose lines are measured on a loaded observation set ``local_continuum``), and
:class:`~.results.ResultReferences` for the document-side half of staleness. What is
genuinely different between the three — their forms, how many steps they have, what a run
synthesises and inverts, their table columns and what they draw on the spectrum — stays in
the controller, through the hooks marked "subclass" below.

Subclasses set the class attributes below (``prefix`` namespaces the DOM/bridge ids,
``kind`` is the left-sidebar focus kind the workspace coordinates, ``ctx.workspace.focus_kind``).
"""
from __future__ import annotations

import html
import time
from copy import deepcopy
from pathlib import Path
from dataclasses import replace

import astropy.units as u
import numpy as np
from bokeh.models import TextInput

from tumult import abundance as ab
from tumult.measured import progress_fraction

from ..execution.jobs import ResultJobs, update_kept_name
from ..execution.processes import format_duration
from ..ui.htmlbridge import pulse, toggle
from ..ui.tableui import SortState, check_header, sort_header
from ..ui.tool import Tool
from ..ui.units import to_range
from ..workspace import store
from .diagnostics import DiagnosticPanels
from .localcontinuum import LocalContinuum
from .results import ResultCollection, ResultReferences

DEFAULT_JUMP_SPAN_AA = 1.0   # jump span for a no-window line when no line has been shown yet
FIT_WINDOW_MARGIN = 0.05     # fraction of a fitted line window shown beyond either edge

# The shared client->server wiring of the three analysis result views, emitted once through
# the ``analysis_script`` template variable. Each tool's markup scopes itself with
# ``data-tool="<prefix>"`` (the step tabs, the results rows, the lines table); the binder
# reads the prefix from the closest scope and writes ``"<kind>:<arg>:<nonce>"`` to the
# ``<prefix>-actions`` channel, which :meth:`ResultView._on_action` splits on colons (uids
# and line indices contain none). A row click jumps to its line; a click on the row's include
# checkbox only toggles the line; a checkbox column's header buttons select all/none of the
# shown rows or filter them (``check-<action>:<column>``). The tools' own run/refit buttons
# and the atom sidebar (which removes its nested corrections) stay in the page script.
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
    var check = ev.target.closest('.check-head-btn[data-checkcol]');
    if (check) {
      analysisAction(prefix, 'check-' + check.dataset.checkact, check.dataset.checkcol);
      return;
    }
    var remaining = ev.target.closest('[data-complete-rungs]');
    if (remaining) { analysisAction(prefix, 'complete-rungs'); return; }
    if (ev.target.closest('[data-update-result]')) { analysisAction(prefix, 'update-result'); return; }
    var reset = ev.target.closest('[data-reset-window]');
    if (reset) { analysisAction(prefix, 'reset-window', reset.dataset.resetWindow); return; }
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
    var cont = ev.target.closest('.line-continuum-chk');
    if (cont) { analysisAction(prefix, 'continuum', cont.dataset.line); return; }
    var chk = ev.target.closest('.line-include-chk');
    if (chk) analysisAction(prefix, 'toggle-line', chk.dataset.line);
  });
})();
"""


def format_value(value, err=None, plus=False, digits=3):
    """Value ± err, or an em-dash when the stage that fills it hasn't run yet."""
    if value is None or not np.isfinite(value):
        return "—"
    sign = "+" if plus else ""
    out = f"{value:{sign}.{digits}f}"
    if err is not None and np.isfinite(err):
        out += f" ± {err:.{digits}f}"
    return out


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


def progress_text(progress, now=None):
    """The strip's right-hand text: steps done (and the running synthesis' own fraction
    when one is under way), time elapsed and the time the remaining work is expected
    to take at the pace of the fraction shown (:func:`tumult.measured.progress_fraction`,
    which only ever advances; nothing to say about that before anything is done).
    The page recomputes it every second from the same attributes (analysis.js), so it
    stays live between snapshots."""
    fraction = progress_fraction(progress)
    elapsed = (time.time() if now is None else now) - float(progress.get("started") or 0)
    text = f"{int(progress.get('done', 0))}/{max(int(progress.get('total', 1)), 1)}"
    if progress.get("run") is not None:
        text += f" · synthesis {100 * float(progress['run']):.0f}%"
    text += f" · {format_duration(elapsed)}"
    if 0 < fraction < 1:
        text += f" · ~{format_duration(elapsed * (1 - fraction) / fraction)} left"
    return text


def progress_html(progress):
    """The strip above a running result's table: the stage under way, a bar of the steps
    done over those expected, and :func:`progress_text`."""
    fraction = progress_fraction(progress)
    pct = 100.0 * fraction
    run = progress.get("run")
    return (f'<div class="result-progress" data-started="{float(progress.get("started") or 0):.3f}" '
            f'data-done="{int(progress.get("done", 0))}" '
            f'data-total="{max(int(progress.get("total", 1)), 1)}" data-fraction="{fraction:.4f}"'
            + (f' data-run="{float(run):.2f}"' if run is not None else "") + '>'
            f'<span class="result-progress-stage">{html.escape(str(progress.get("stage", "")))}</span>'
            f'<div class="result-progress-bar"><div class="result-progress-fill" '
            f'style="width: {pct:.1f}%"></div></div>'
            f'<span class="result-progress-text">{html.escape(progress_text(progress))}</span></div>')


class ResultView(Tool):
    # --- subclass configuration --------------------------------------------------------
    prefix = ""            # DOM/bridge namespace: "abund" | "corr" | "ionbal"
    steps = 3              # how many step panels the view's switcher drives
    kind = ""              # the workspace focus_kind value while a result owns the centre view
    section = ""           # store section for the view's own form values
    empty_table_msg = "No result selected"
    y_axis_label = "A(X)"
    local_continuum = False   # subclass: True where a line carries a continuum of its own

    # sortable per-line columns shared by every result table: {col: sort-value(ln)}. The
    # header markup lives in each subclass's table_html; here we own the state + the row
    # order so the tables sort identically (and stay in step with the Line Navigator).
    # Wavelength / excitation read low→high, EW and abundance magnitude are most useful
    # largest-first.
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
        self.setup = setup                # SynthSetup: the view's own synthesis config
        self.own_mode_bridges((".sidebar-right", f"{self.prefix}-active"))
        self._view_bridge = ctx.bridges.cls("#main", f"{self.prefix}-mode")
        self._step_bridges = [ctx.bridges.cls("#main", f"{self.prefix}-step-{i}")
                              for i in range(1, self.steps + 1)]
        self._results_bridge = ctx.bridges.html(f"{self.prefix}-results")
        self._table_bridge = ctx.bridges.html(f"{self.prefix}-table")
        self._progress_bridge = ctx.bridges.html(f"{self.prefix}-progress")
        # another list took (or gave back) the single left-sidebar selection: restyle our rows
        ctx.subscribe("focus-changed", lambda **_: self.push_results())

        self.step = 1
        self.check_filters = set()         # checkbox columns showing only checked lines
        self.results = ResultCollection()
        self.references = ResultReferences(self.canvas, self._result_spec)
        self.jobs = ResultJobs(ctx, self.find, self._job_progress, self._job_run_progress)
        if self.task_kind:
            # a record of this kind another document submitted lands here on a reconnect
            ctx.task_importers[self.task_kind] = self.import_task_result
            ctx.task_progress[self.task_kind] = self.import_task_progress
        self._last_span_aa = None          # span (Å) of the last shown line, for no-window jumps
        self.jumped = None                 # index of the line whose window is accented
        self._opened_on_row = False        # the selected result was shown on its first row

        # step switch / run buttons / result clicks / line toggles all report through one
        # client->server channel, written as ``"<kind>:<arg>:<nonce>"``.
        self.actions = TextInput(value="", visible=False, name=f"{self.prefix}-actions")
        self.actions.on_change("value", self._on_action)
        self.continuum = LocalContinuum(
            ctx, self.prefix, spec_of=self._result_spec, draw_hatch=self._set_hatch,
            on_edited=self._continuum_edited, current=lambda: (self.selected(), self.jumped),
            listed=lambda res: self.find(res.uid) is res) if self.local_continuum else None

    # --- the diagnostic panels ------------------------------------------------------------
    def build_panels(self, source=None):
        """Create this view's pair of diagnostic panels (subclasses call it once they have
        their own sources, then add whiskers and any extra series)."""
        self.panels = DiagnosticPanels(self.prefix, self.y_axis_label, source)
        return self.panels

    @property
    def figs(self):
        return self.panels.figs

    @property
    def src(self):
        return self.panels.src

    # --- window highlights ---------------------------------------------------------------
    # A result's windows belong to the result, not to a jump: a view that paints them shows
    # them for as long as it owns the centre view, with the jumped-to line accented. Only the
    # corrections view lets that accent be dragged (its fit windows are its own to edit);
    # the others' windows are edited in the LineMasks view. The call sites are the
    # subclass's, so a view whose lines carry no user-drawn windows paints nothing.
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
    # Like the windows, the syntheses belong to the result: a view that keeps one per line
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
        return self.ctx.workspace.focus_kind == self.kind

    def push_view(self, state):
        """Show/hide the view's centre pane and hand the selection focus around."""
        pulse(self._view_bridge, state)
        self.ctx.workspace.focus(self.kind if state == "on" else "spectrum")

    def set_step(self, n):
        self.step = max(1, min(self.steps, int(n)))
        for i, bridge in enumerate(self._step_bridges, start=1):
            toggle(bridge, i == self.step)

    # --- Tool protocol -------------------------------------------------------------------
    @property
    def models(self):
        return [*self.form1.models, self.actions, *self.setup.models]

    def activate(self):
        # first open: seed the view's own synthesis setup from the synthesis modal, then
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
        if self.continuum is not None:
            self.continuum.hide()
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
            self.ctx.removal.request(self.name, arg)
        elif kind == "toggle-line":
            self._toggle_line(arg)
        elif kind == "continuum":
            self._toggle_continuum(arg)
        elif kind == "complete-rungs":
            self._complete_rungs()
        elif kind == "update-result":
            self.update_result()
        elif kind == "check-filter":
            if arg in self.check_cols:
                self.check_filters ^= {arg}
                self._render_selected()
        elif kind in ("check-all", "check-none"):
            self._set_checks(arg, kind == "check-all")
        elif kind == "jump-line":
            self.jump_line(arg)
        elif kind == "sort":
            if self.sort.toggle(arg):
                res = self.selected()
                if res is not None:
                    self.push_table(self.table_html(res))
        else:
            self._extra_action(kind, arg)

    def _extra_action(self, kind, arg):
        """Subclass hook: the view-specific run actions."""

    # --- staleness: a solution whose inputs have changed since ------------------------------
    # A result is out of date when what it was solved from differs from what it would be
    # solved from now. The scientific half of that comparison is the result type's own
    # (``stale_reasons`` of the result); the document half — which spectrum its reference
    # names, its digest now, and how the recorded preparation differs from the current one —
    # is :class:`~.results.ResultReferences`. Nothing is recomputed; the panel says what
    # changed and offers ``update_label``.
    update_label = ""      # subclass: the action that solves the selected result again

    def update_result(self):
        """Subclass hook: solve the selected result again from its current inputs."""

    def stale_reasons(self, res):
        """``{key: text}`` saying why the solved ``res`` no longer matches its inputs; empty
        before it is solved and while it is being computed."""
        if res is None or res.stage != self.final_stage or self.is_computing(res):
            return {}
        return self.references.reasons(res)

    def sync_target(self, res, spec):
        """Before ``res`` is solved again on ``spec``: when its target changed, measure the
        lines' EWs again on ``spec`` and record the current preparation, so the new solution
        stands on the target as it is now."""
        if "target" not in self.stale_reasons(res):
            return
        outcome = ab.refresh_ews(spec, res.lines, profile=res.profile, snr=res.snr,
                                 cont_pct=res.cont_pct)
        line_list = (res.preparation or {}).get("line_list")
        res.preparation = self.canvas.preparation(self.canvas.spectra.find(res.ref_uid))
        if line_list:
            res.preparation["line_list"] = line_list
        self.report_refresh(res, outcome)

    def stale_banner(self, res):
        """The panel's statement that ``res`` shows its previous solution: what changed and
        the actions that bring it up to date."""
        reasons = self.stale_reasons(res)
        if not reasons:
            return ""
        items = "".join(f"<li>{html.escape(text)}</li>" for text in reasons.values())
        rungs = ('<button class="text-btn" data-complete-rungs>Compute remaining abundance '
                 'rungs</button>' if any(ln.missing_rungs for ln in res.lines) else "")
        label = html.escape(self.update_label)
        return (f'<div class="stale-banner" data-tool="{self.prefix}" role="status">'
                '<div class="stale-banner-head"><i class="ti ti-alert-triangle"></i> '
                'Out of date: the values below are the previous solution</div>'
                f'<ul>{items}</ul>'
                f'<button class="primary-btn" data-update-result>{label}</button>{rungs}'
                f'<div class="hint">{label} replaces this solution; use Edit → Duplicate '
                'first to keep it.</div></div>')

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
    def find(self, uid):
        return self.results.find(uid)

    def selected(self):
        return self.results.selected_result()

    # --- durable calculations: this view's task kind -----------------------------------
    # The view's calculation runs as a workflow record of ``task_kind``
    # (tumult.synthesis.tasks), so it outlives the tab and the server. Its snapshots
    # and result are result objects carrying the uid of the one they continue; a uid this
    # list does not hold (a record another document submitted) is adopted as a new
    # result. Subclasses that submit through ``self.jobs.submit_operation`` name the operation.
    task_kind = None

    def _adopt(self, incoming, name):
        """The listed result ``incoming`` continues, updated in place from it (its name
        stays), or ``incoming`` itself, listed as ``name`` and selected."""
        res = self.find(incoming.uid)
        if res is None:
            self._add_result(incoming, name)
            self.push_results()
            return incoming
        update_kept_name(res, incoming)
        return res

    def _add_result(self, res, name):
        """List ``res`` under ``name``, made unique in the document, and select it."""
        res.name = self.ctx.workspace.unique_name(name, pending=False)
        self.results.append(res)

    def import_task_progress(self, snapshot, handle):
        self._job_progress(self._adopt(snapshot, handle.label))

    def import_task_result(self, result, handle):
        self._job_progress(self._adopt(result, handle.label))

    # --- session persistence: the envelope every result list shares ----------------------
    # The result objects themselves are built and relinked by :mod:`.results`, without a
    # view; only the list/selected-uid envelope and what it takes to show it live here.
    def snapshot(self):
        """``{"results": [...], "selected": uid}`` for the session slot."""
        return self.results.snapshot()

    def restore(self, results, selected):
        """Replace the result list (and repaint the sidebar); ``selected`` must be one of
        them or ``None``. Pending result jobs of the old list are forgotten."""
        self.jobs.current.clear()
        for res in results:
            res.progress = {}            # nothing of this document works on it (yet)
        self.results.replace(results, selected)
        self.jumped = None
        self._opened_on_row = False
        self.push_results()

    def select_result(self, uid):
        if self.find(uid) is None:
            return
        self.results.selected = uid
        self.ctx.tools.activate(self.name)   # a result click opens the panel (no toolbar needed)
        self._opened_on_row = False
        self._clear_line_decor()
        self.push_results()
        self._render_selected()
        self.push_view("on")

    def blur(self):
        """A spectrum object took focus: leave the centre view but keep the selection."""
        self.push_view("off")
        self._clear_line_decor()
        self.push_results()                  # drop the result highlight

    def is_computing(self, res):
        """Whether a submission of this document still works on ``res``: one bound to it
        (keyed by its uid) or the one that creates it and delivers its later stages
        (keyed by the name it was submitted under)."""
        return res.uid in self.jobs.current or res.name in self.jobs.current

    def remove_result(self, uid):
        """Drop result ``uid`` and cancel the workflow computing it. Returns what
        :meth:`insert_result` needs to put it back, or ``None``."""
        res = self.find(uid)
        if res is None:
            return None
        placed = {"res": res, "index": self.results.index(res),
                  "selected": self.results.selected == uid}
        self.jobs.invalidate(uid)
        self.jobs.invalidate(res.name)
        was_selected = self.results.remove(res)
        self.push_results()
        if self.results.selected and was_selected:
            self.select_result(self.results.selected)
        elif not len(self.results):
            self.push_view("off")
            self._clear_line_decor()
        return placed

    def insert_result(self, res, index, selected):
        """Put back a result :meth:`remove_result` took out, at its position and with its
        uid and state; ``selected`` shows it again. Its cancelled workflow stays cancelled."""
        self.results.insert(res, index)
        if selected:
            self.select_result(res.uid)
        else:
            self.push_results()

    def duplicate_result(self, uid=None):
        """Add a copy of result ``uid`` (the selected one by default) under a unique name
        and select it. The copy has its own uid, so the original's workflow keeps reaching
        the original only, and a rerun or curation of the copy leaves the original as it
        is. A result still being computed is not copied."""
        res = self.find(uid or self.results.selected)
        if res is None:
            return None
        if self.is_computing(res):
            self.notify(f"{res.name} is still being computed; duplicate it once it has "
                        "finished or been cancelled.")
            return None
        copy = self.results.duplicate(res)
        self._add_result(copy, res.name)
        self.select_result(copy.uid)
        self.notify(f"Duplicated as “{copy.name}”.", level="success")
        return copy

    def rename_result(self, uid, name):
        """Give result ``uid`` the name ``name``, made unique in the document. A result still
        being computed keeps its name: its workflow delivers to the name it was submitted
        under."""
        res = self.find(uid)
        if res is None:
            return None
        if self.is_computing(res):
            self.notify(f"{res.name} is still being computed; rename it once it has "
                        "finished or been cancelled.")
            return None
        res.name = self.ctx.workspace.unique_name(name, own=res.name)
        self.push_results()
        return res

    def _toggle_line(self, token):
        res = self.selected()
        if res is None:
            return
        try:
            idx = int(token)
        except ValueError:
            return
        if 0 <= idx < len(res.lines):
            self._set_include(res, [res.lines[idx]], not res.lines[idx].include)

    def _set_include(self, res, lines, on):
        self.jobs.invalidate(res.uid)
        for ln in lines:
            ln.include = on
        if res.stage == self.final_stage:
            res.summary = self._summarise(res.lines)
            self.push_results()
        self._render_selected()

    def _set_checks(self, col, on):
        """A header's select all/none: set checkbox column ``col`` on every shown line."""
        res = self.selected()
        if res is None or col not in self.check_cols:
            return
        lines = [res.lines[i] for i in self._line_order(res)
                 if self.CHECK_VALUES[col](res.lines[i]) != on]
        if not lines:
            return
        if col == "include":
            self._set_include(res, lines, on)
        elif self.continuum is not None:
            self.continuum.set_enabled(res, lines, on)

    # subclass: the stage name after the last run, the aggregation and the sidebar summary
    final_stage = "uncert"

    def _summarise(self, lines):
        raise NotImplementedError

    def _result_summary(self, res):
        """The sidebar's second line: the result's current value or stage."""
        raise NotImplementedError

    # --- streaming + rendering -------------------------------------------------------------
    def _job_progress(self, res):
        self.push_results()
        self._stream_row(res)

    def _job_run_progress(self, res):
        """Doc-thread callback every few seconds while ``res``'s dispatch run advances:
        only its progress indication is repainted, the table and panels are not."""
        self.push_results()
        if self.selected() is res and self.is_focused():
            self.push_progress(res)

    def _stream_row(self, res):
        """Doc-thread callback: repaint everything an incoming line changed — the table, the
        diagnostic panels and the synthetic overlays — so a fit is watched as it runs."""
        if self.selected() is res and self.is_focused():
            self._ensure_jumped(res)
            self.push_table(self.table_html(res))
            self._render_scatter()
            self._render_overlays(res)

    def _ensure_jumped(self, res):
        """Show a result on its first row: the view's line panels (the curve of growth,
        the hatch, the accented window) all belong to the jumped-to line, so a result
        whose lines have arrived opens on the first row of the table as displayed rather
        than on nothing. Done once per selection; clearing the decor afterwards (a
        no-window line, say) leaves no row selected until the next click."""
        if self.jumped is not None or self._opened_on_row or not res.lines:
            return
        self._opened_on_row = True
        self.jump_line(self._line_order(res)[0])

    def submit(self, res, *args, **kwargs):
        """Submit ``res``'s calculation (:meth:`ResultJobs.submit_operation`, with
        ``result=res``) and repaint it as a running result: the sidebar row and the
        table pane show its progress strip from here on."""
        job = self.jobs.submit_operation(*args, result=res, **kwargs)
        self._job_progress(res)
        return job

    def progress_strip(self, res):
        """The progress strip for ``res``, or nothing once its record is cleared: by the
        operation with its result, by the job when it stops, by a session load."""
        progress = getattr(res, "progress", None) or {}
        return progress_html(progress) if progress else ""

    def push_results(self):
        self._results_bridge.value = self.results_html()

    def push_progress(self, res):
        """Show the selected result's progress strip over the table, or nothing."""
        self._progress_bridge.value = self.progress_strip(res) if res is not None else ""

    def push_table(self, html_text):
        """Show ``html_text`` in the bottom pane, with the selected result's progress
        strip over it while it is being computed."""
        self.push_progress(self.selected())
        self._table_bridge.value = html_text

    # --- shared sortable-table scaffolding ------------------------------------------------
    # the checkbox columns a table shows, and each one's per-line state
    check_cols = ("include",)
    CHECK_VALUES = {"include": lambda ln: ln.include, "continuum": lambda ln: ln.renorm}

    def _line_order(self, res):
        """Indices into ``res.lines`` of the shown rows, in the active sort order (NaN/None
        keys sort last). A filtered checkbox column hides the lines it leaves unchecked.

        Rows carry their *original* index as ``data-line`` so the client's toggle/jump
        actions keep addressing ``res.lines[i]`` no matter how the display is sorted.
        """
        order = self.sort.order(res.lines, self.sort_keys[self.sort.col])
        filters = [self.CHECK_VALUES[col] for col in self.check_filters]
        return [i for i in order if all(f(res.lines[i]) for f in filters)]

    def _table_head(self, res, columns):
        """The clickable ``<thead>`` row: the leading include column, then ``columns``.

        ``columns`` is a sequence of ``(label, col)`` pairs. A ``col`` in :attr:`check_cols`
        renders a checkbox-column header, a ``None`` column (e.g. the flags column) a plain,
        non-sortable one, and any other a sortable one.
        """
        shown = [res.lines[i] for i in self._line_order(res)]
        cells = []
        for label, col in (("", "include"), *columns):
            if col in self.check_cols:
                cells.append(check_header(label, col, [self.CHECK_VALUES[col](ln) for ln in shown],
                                          col in self.check_filters, noun="lines"))
            else:
                cells.append(sort_header(label, col, self.sort))
        return "<tr>" + "".join(cells) + "</tr>"

    def _filtered_out_row(self, columns):
        """A body row saying the column filters hide every line; the header stays so the
        filters can be switched off again."""
        if not self.check_filters:
            return ""
        return (f'<tr><td class="linelist-empty" colspan="{len(columns) + 1}">'
                'No line passes the column filters.</td></tr>')

    def _render_selected(self):
        res = self.selected()
        if res is None:
            self._render_scatter()
            self.push_table(f'<div class="linelist-empty">{self.empty_table_msg}</div>')
            return
        self._ensure_jumped(res)
        self._render_scatter()
        self.push_table(self.table_html(res))

    # --- jump-to-line -----------------------------------------------------------------------
    def jump_line(self, token):
        """Show line ``token`` of the selected result: move the view onto its window and
        paint the decor that belongs to it."""
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

    # --- the left sidebar's rows ---------------------------------------------------------
    def _results_rows(self, badge, results=None, nested=False):
        """Sidebar rows for ``results`` (default: all); ``nested`` indents them under the
        object they were computed from. The row carries ``data-tool`` (this view's
        prefix) and ``data-result`` (the uid), its button ``data-result-remove``: the
        attributes the page's analysis binder routes select/remove clicks by."""
        focused = self.is_focused()
        rows = []
        nest = " result-nested" if nested else ""
        for r in (self.results if results is None else results):
            sel = " selected" if (r.uid == self.results.selected and focused) else ""
            summary = self._result_summary(r)
            progress = getattr(r, "progress", None) or {}
            if progress:
                summary = (f"{progress.get('stage', 'running')} · "
                           f"{100 * progress_fraction(progress):.0f}%")
            rows.append(
                f'<div class="object-item result-item{nest}{sel}" '
                f'data-tool="{self.prefix}" data-result="{r.uid}">'
                f'<span class="result-badge">{badge}</span>'
                f'<span class="result-meta"><span class="object-name" '
                f'title="{html.escape(r.name)}">{html.escape(r.name)}</span>'
                f'<span class="result-summary">{html.escape(summary)}</span></span>'
                + ('<span class="lm-dirty" title="Its inputs changed after it was solved">'
                   'Out of date</span>' if self.stale_reasons(r) else "") +
                f'<button class="result-btn" data-result-remove="{r.uid}" '
                f'title="Remove from workspace"><i class="ti ti-x"></i></button></div>')
        return rows

    # --- measurement-owned local continua ------------------------------------------
    def _set_hatch(self, res, ln):
        """Hatch the area ``ln``'s measured EW stands for, through its local continuum when
        enabled. A line restored from a task or session has no live profile; it is measured
        again here.

        The hatch spans :func:`tumult.abundance.ew_window`, not the fit window: a fitted
        profile's width counts wings that run past the window it was fitted over, so it is
        drawn right across the view and follows it (``canvas.set_abund_hatch``). Only a
        directly integrated width stops at the window, because that is all it integrated.
        """
        if ln.fit_func is None and np.all(np.isfinite(ln.fit_window)):
            spec = self._result_spec(res)
            if spec is not None:
                try:
                    ln.fit_func = ab.measured_profile_func(ln, spec, res.profile)
                except (ValueError, RuntimeError):
                    pass
        self.canvas.set_abund_hatch(ln.fit_func, ab.ew_window(ln, res.profile),
                                    nodes=ln.continuum if ln.renorm else None)

    def _render_continuum(self):
        """Show the local-continuum markers of the line in view, or hide them."""
        if self.continuum is None:
            return
        res = self.selected()
        if (not self.active or not self.is_focused() or res is None or self.jumped is None
                or self.jumped >= len(res.lines)):
            return self.continuum.hide()
        self.continuum.render(res, self.jumped)

    def _continuum_edited(self, res, indices):
        """A local-continuum measurement landed on lines ``indices`` of ``res``: the
        solution standing on them is superseded, so the pending run is dropped and
        everything that shows them is repainted."""
        self.jobs.invalidate(res.uid)
        self.push_results()
        self._render_selected()
        if self.jumped is not None and self.jumped in indices:
            self._jump_decor(res, res.lines[self.jumped])
        self.ctx.workspace_changed()

    def _toggle_continuum(self, token):
        res = self.selected()
        if res is None or self.continuum is None:
            return
        idx = int(token)
        if not 0 <= idx < len(res.lines):
            return
        self.jumped = idx
        self.continuum.edit(res, idx, enabled=not res.lines[idx].renorm)
        self.jump_line(token)

    def _complete_rungs(self):
        res = self.selected()
        if res is None or not any(ln.missing_rungs for ln in res.lines):
            return
        target = deepcopy(self._result_spec(res))

        def work(working, token, progress):
            return ab.complete_result_rungs(working, target=target,
                                            cancel_token=token, progress=progress)

        def done(result):
            self.push_results()
            self._render_selected()
            self.notify("Cached curves of growth updated." if not any(
                ln.missing_rungs for ln in result.lines) else
                "More abundance rungs are still needed; inspect the flagged lines.")

        self.jobs.submit(f"rungs:{res.name}", work, done, result=res, detached=True)
