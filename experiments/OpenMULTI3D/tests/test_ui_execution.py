"""Background computations must not apply stale UI state or wait on HOP locks."""
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from bokeh.models import ColumnDataSource, Range1d, TapTool

from tumult import Spectrum
from tumult.synthesis.cancel import CancelToken


@pytest.fixture(autouse=True)
def _temporary_config(tmp_path, monkeypatch):
    """Every document built here reads and writes a config of its own; the project's
    ``.tumultstudio.json`` is neither consulted nor polluted."""
    import store
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "tumultstudio.json")


def test_scatter_ranges_include_error_bars():
    from abundwizard import AbundanceWizard, PANEL_AXES

    wizard = AbundanceWizard.__new__(AbundanceWizard)
    wizard.figs = []
    wizard.mean_lines = []
    wizard.mean_bands = []
    for key in PANEL_AXES:
        wizard._new_fig(key)
    points = {
        "x_evi": [1.0, 2.0],
        "x_rew": [-5.5, -4.5],
        "a": [0.0, 1.0],
        "lo": [-3.0, 0.5],
        "hi": [0.25, 4.0],
    }

    wizard._crop_panels((points,))

    for panel in wizard.figs:
        assert isinstance(panel.y_range, Range1d)
        assert panel.y_range.start == pytest.approx(-3.35)
        assert panel.y_range.end == pytest.approx(4.35)
        assert panel.y_range.reset_start == pytest.approx(-3.35)
        assert panel.y_range.reset_end == pytest.approx(4.35)


def test_abundance_line_jump_pads_fit_window():
    from abundwizard import AbundanceWizard

    line = SimpleNamespace(window=(5000.0, 5002.0))
    result = SimpleNamespace(lines=[line])
    ranges = []
    decorated = []
    wizard = AbundanceWizard.__new__(AbundanceWizard)
    wizard._last_span_aa = None
    wizard.selected = lambda: result
    wizard._set_range = lambda lo, hi: ranges.append((lo, hi))
    wizard._jump_decor = lambda res, ln: decorated.append((res, ln))

    wizard._jump_line("0")

    assert ranges == [pytest.approx((4999.9, 5002.1))]
    assert decorated == [(result, line)]


def test_scatter_click_selects_and_reveals_table_row():
    from abundwizard import AbundanceWizard, PANEL_AXES

    wizard = AbundanceWizard.__new__(AbundanceWizard)
    wizard.prefix = "abund"
    wizard.figs = []
    wizard.mean_lines = []
    wizard.mean_bands = []
    for key in PANEL_AXES:
        wizard._new_fig(key)
    source = ColumnDataSource({f"x_{key}": [] for key in PANEL_AXES} | {
        "a": [], "idx": [],
    })

    wizard._add_series(source, "#ffffff")

    callback = source.selected.js_property_callbacks["change:indices"][0]
    assert "table.scrollTop" in callback.code
    assert "row.click()" in callback.code
    assert all(any(isinstance(tool, TapTool) for tool in panel.tools)
               for panel in wizard.figs)
    for panel in wizard.figs:
        renderer = panel.renderers[-1]
        assert renderer.nonselection_glyph.fill_alpha == 1.0
        assert renderer.nonselection_glyph.line_alpha == 1.0


def test_band_ready_does_not_wait_for_native_calculation():
    from tumult import linestrength
    from tumult.synthesis.runtime import HOP_LOCK
    star = SimpleNamespace(key='test-pending')
    with HOP_LOCK, ThreadPoolExecutor(1) as worker:
        assert worker.submit(linestrength.band_ready, star, 'none', 1, 2).result(1) is False


def test_annotation_folder_omits_molecules_only_above_view_limit(tmp_path):
    from pathlib import Path
    from tumult.linestrength import LineStrengths, MOLECULE_VIEW_LIMIT_AA

    def line_file(name, species):
        path = tmp_path / name
        path.write_text(f"'{species}'\n'0' 0\n")
        return path

    atom = line_file("atoms.txt", "Fe I")
    hydrogen = line_file("hydrogen", "H I  descriptive header")
    line_file("tio.txt", "TiO")
    strengths = LineStrengths(master_path=atom, linelist=tmp_path, refstar="gdwarf")
    obj = SimpleNamespace(meta={})

    _star, narrow, _band = strengths.request_for(
        obj, 5000.0, 5000.0 + MOLECULE_VIEW_LIMIT_AA)
    assert Path(narrow) == tmp_path

    _star, wide, _band = strengths.request_for(
        obj, 5000.0, 5000.0 + MOLECULE_VIEW_LIMIT_AA + 0.01)
    assert tuple(path.name for path in wide) == (atom.name, hydrogen.name)

    catalog = tmp_path.parent / "catalog.txt"
    catalog.write_text(
        "'Fe I'\n'0' 1\n5000.0 1.0 -1.0 0 0 0 0 0\n"
        "'TiO'\n'0' 1\n5001.0 1.0 -1.0 0 0 0 0 0\n")
    strengths.master_path = str(catalog)
    source = strengths.lines_for(obj, 5000.0, 5000.0 + MOLECULE_VIEW_LIMIT_AA + 0.01)
    assert source.lines.ion.tolist() == ["Fe I"]


def test_normalisation_uses_snapshot_and_discards_stale_completion():
    from normalisetool import NormaliseTool
    tool = NormaliseTool.__new__(NormaliseTool)
    obj = SimpleNamespace(raw=Spectrum(np.arange(5000., 5010., .01), np.ones(1000)),
                          already_normalised=False)
    values = dict(res=0., win=0., pct=0., peak=1.)
    queued, rendered = [], []
    tool.active = True
    tool._selection_request = 0
    tool._obs = lambda: obj
    tool.form = SimpleNamespace(value=lambda k: values[k])
    tool.knots = SimpleNamespace(data={'x': [], 'y': []})
    tool.canvas = SimpleNamespace(set_norm_view=rendered.append)
    tool.ctx = SimpleNamespace(submit=lambda name, run, done, detached=False: queued.append((run, done)))

    tool._render_selection()
    assert not rendered  # submitting work does not change the view
    run, done = queued.pop()
    values['pct'] = 50.  # simulate an edit while the computation is queued
    with ThreadPoolExecutor(1) as worker:
        result = worker.submit(run, CancelToken()).result(3)
    done(result)
    assert not rendered

    tool._render_selection()
    run, done = queued.pop()
    with ThreadPoolExecutor(1) as worker:
        result = worker.submit(run, CancelToken()).result(3)
    assert not rendered  # only the document-thread completion changes the view
    done(result)
    assert len(rendered) == 1


def test_ion_balance_setup_is_private_and_interpolation_only(monkeypatch):
    import ionbalance

    grids = SimpleNamespace(
        names=[("marcs", "MARCS"), ("stagger", "STAGGER")],
        _find=lambda gid: ({"id": gid, "name": gid.upper(), "folder": f"/{gid}",
                            "pattern": "pattern", "format": "Marcs"}
                           if gid in {"marcs", "stagger"} else None),
    )
    # Even a globally selected individual atmosphere must have no influence on this tool.
    sections = {
        "ionbal": {"atmos_grid_id": "stagger"},
        "synthesis": {"atmos_grid_id": "__individual__", "ref_grid_id": "marcs"},
        "synth_run": {"atmos_mode": "pick", "atmos_file": "/fixed.mod",
                      "mode": "NLTE", "abundances": "C=9"},
    }
    monkeypatch.setattr(ionbalance.store, "get_section", lambda name: dict(sections.get(name, {})))
    monkeypatch.setattr(ionbalance, "load_prefs", lambda: {"linelist": "lines", "abund_file": "abund"})
    monkeypatch.setattr(ionbalance, "missing_prefs", lambda prefs: [])

    setup = ionbalance.BalanceSynthesis(grids)
    context = setup.context()

    assert setup.form.value("atmos_grid_id") == "stagger"
    assert [value for value, _ in setup.form.fields["atmos_grid_id"].options] == [
        "marcs", "stagger"]
    # a SynthRequest, like synthsetup.context(): solve_ionisation_balance reads attributes,
    # and the dict this used to return made the solve fail on request.prefs
    assert context.grid["id"] == "stagger"
    assert context.atmos_file is None
    assert context.run_vals["atmos_mode"] == "interp"
    assert context.run_vals["mode"] == "LTE"
    assert context.run_vals["abundances"] == ""


def test_synthesis_geometry_is_inferred_from_model_format(tmp_path):
    from synthsetup import model_is_3d

    assert not model_is_3d("Marcs")
    assert not model_is_3d("Text")
    assert model_is_3d("Stagger")
    assert model_is_3d("Co5bold")
    assert model_is_3d("M3DIS")

    saved = tmp_path / "saved"
    saved.mkdir()
    (saved / "atmos_mesh.txt").write_text("1 1 128\n")
    assert not model_is_3d("save", saved)
    (saved / "atmos_mesh.txt").write_text("8 8 128\n")
    assert model_is_3d("save", saved)


def test_result_jobs_snapshot_progress_and_reject_stale_completion():
    from jobs import ResultJobs
    queued, callbacks, rendered, completed = [], [], [], []
    result = SimpleNamespace(uid='fit', lines=[SimpleNamespace(a=0.)], termination='measured')
    results = {'fit': result}
    ctx = SimpleNamespace(submit=lambda name, run, done, detached=False: queued.append((run, done)),
                          schedule=callbacks.append, notify=lambda *a: None, closed=False)
    jobs = ResultJobs(ctx, results.get, lambda r: rendered.append(r.lines[0].a))

    def work(working, token, progress):
        working.lines[0].a = 1.
        progress(working)
        working.lines[0].a = 2.
        progress(working)
        return working

    jobs.submit('fit', work, completed.append, result=result)
    run, done = queued.pop()
    outcome = run(CancelToken())
    assert result.lines[0].a == 0.
    assert len(callbacks) == 1  # one pending paint, independent of worker updates
    callbacks.pop()()
    assert result.lines[0].a == 1.
    done(outcome)
    assert result.lines[0].a == 2. and completed == [result]

    jobs.submit('fit', work, completed.append, result=result)
    run, done = queued.pop()
    outcome = run(CancelToken())
    results.pop('fit')  # removal while running cannot resurrect the object
    done(outcome)
    for callback in callbacks:
        callback()
    assert completed == [result]


def test_result_jobs_isolate_errors_cancel_and_closed_documents():
    from jobs import ResultJobs
    queued, callbacks, notes, completed = [], [], [], []
    result = SimpleNamespace(uid='fit', termination='measured', converged=True)
    ctx = SimpleNamespace(submit=lambda name, run, done, detached=False: queued.append((run, done)),
                          schedule=callbacks.append, notify=notes.append, closed=False)
    jobs = ResultJobs(ctx, lambda uid: result, lambda r: None)

    def fail(working, token, progress):
        raise ValueError('this job failed')

    jobs.submit('first', fail, completed.append, result=result)
    run, done = queued.pop()
    done(run(CancelToken()))
    assert notes == ['first failed: this job failed']
    assert result.termination == 'failed' and not result.converged

    def cancel(working, token, progress):
        token.cancel()
        token.check()

    jobs.submit('second', cancel, completed.append, result=result)
    run, done = queued.pop()
    from tumult.synthesis.cancel import CancelledError
    with pytest.raises(CancelledError):
        run(CancelToken())
    callbacks.pop()()
    assert result.termination == 'cancelled' and not completed

    jobs.submit('third', lambda working, token, progress: working, completed.append,
                result=result)
    run, done = queued.pop()
    outcome = run(CancelToken())
    ctx.closed = True
    done(outcome)
    assert not completed


def test_restored_abundance_uses_identity_with_duplicate_names():
    from abundancetool import AbundanceTool
    from tumult.abundance import AbundanceResult
    tool = AbundanceTool.__new__(AbundanceTool)
    tool.list_observations = lambda: [('one', 'same'), ('two', 'same')]
    tool.push_results = lambda: None
    result = AbundanceResult(element='Fe', symbol='Fe', summary={}, ref_uid='two',
                             ref_name='same', uid='result')
    tool.restore({'results': [result.to_dict()]})
    assert tool.results[0].ref_uid == 'two'
    legacy = result.to_dict()
    legacy.pop('schema_version')
    legacy['ref_uid'] = 'old'
    tool.restore({'results': [legacy]})
    assert tool.results[0].ref_uid is None  # ambiguous names do not invent an association


def _build_document():
    """Build the real TumultStudio2 document headlessly; return it with its context."""
    from bokeh.document import Document
    from bokeh.io.state import curstate

    document = Document()
    curstate().document = document
    import app
    return document, app.build()


def test_closing_a_session_spares_only_detached_jobs():
    """A tab closing cancels this document's work, except the jobs that asked to survive:
    a multi-hour departure solve cannot assume the tab stays open."""
    import threading

    document, ctx = _build_document()
    tokens, started, release = {}, threading.Event(), threading.Event()

    def loader(name):
        def work(token):
            tokens[name] = token
            if len(tokens) == 2:
                started.set()
            release.wait(5)
        return work

    ctx.submit("keeper", loader("keeper"), lambda *_: None, detached=True)
    ctx.submit("victim", loader("victim"), lambda *_: None)
    assert started.wait(5), "both jobs must be running before the session closes"

    for callback in list(document.session_destroyed_callbacks):
        callback(None)          # what a closing browser tab fires

    assert not tokens["keeper"].cancelled
    assert tokens["victim"].cancelled
    release.set()


def test_mask_tool_range_feeds_line_masks_and_the_synthesis_modal():
    """One drawn range starts a line mask, extends it through "+ New mask", and opens the
    synthesis modal on either the free range or the set's windows (docs/GUIWorkflow.md)."""
    document, ctx = _build_document()
    messages = []
    ctx.notify = lambda text, level="info": messages.append(text)
    x = np.arange(5000., 5010., .01)
    ctx.canvas.add_object(Spectrum(x, 1 - 0.5 * np.exp(-0.5 * ((x - 5005) / 0.1) ** 2)),
                          name="obs", source=("file", "/nowhere/obs.fits"))
    tools = ctx.tools
    range_tool, masks_view = tools._tools["range"], tools._tools["linemask"]

    tools.toggle("range")
    range_tool._on_select(5004.5, 5005.5)
    assert range_tool.selection == (5004.5, 5005.5)
    panel = range_tool.panel_html()
    assert "range-new-mask" in panel and "ew-method-select" not in panel
    range_tool.toggle_ew()                      # the EW is measured only on request
    assert range_tool.ew._weq is not None and "ew-method-select" in range_tool.panel_html()
    range_tool.toggle_ew()                      # ... and hidden again
    assert range_tool.ew._weq is None and "ew-method-select" not in range_tool.panel_html()

    range_tool.new_mask()                       # always a fresh set, opened in its view
    assert tools.active_name() == "linemask"
    assert [(s.name, len(s.masks)) for s in ctx.masks.listed] == [("New line mask", 1)]

    masks_view.actions.value = "new 1"           # "+ New mask" reuses the gesture in add mode
    assert tools.active_name() == "range" and range_tool._add[0] == "New line mask"
    range_tool._on_select(5001.0, 5002.0)
    assert "range-add" in range_tool.panel_html()
    range_tool.add_window()
    assert tools.active_name() == "linemask"
    assert [(m.lmin, m.lmax) for m in ctx.masks.masks] == [(5004.5, 5005.5), (5001.0, 5002.0)]

    modal = masks_view.on_synth.__closure__[0].cell_contents
    modal.open_range(5004.5, 5005.5, "synth")
    assert modal.window_form.values() == {"xmin": 5004.5, "xmax": 5005.5}
    assert "Synthesise" in modal.actions_html()
    modal.open_masks(ctx.masks.mset, "fit")
    assert "Fit" in modal.actions_html()
    # the fit synthesises the span of the synth windows and compares inside the fit windows
    assert modal._fit_windows() == (5000.5, 5006.0, [(5004.5, 5005.5), (5001.0, 5002.0)])


# --- durable synthesis workflows on the document ---------------------------------------
def _stub_request(root):
    """A resolvable synthesis request whose execution the tests stub out."""
    from tumult.synthesis.request import SynthRequest
    return SynthRequest(prefs={"abund_file": None, "dispatch_cwd": str(root)},
                        grid={"folder": ""},
                        run_vals={"abundances": "", "relative_abundances": False,
                                  "mode": "LTE"},
                        target=(5750.0, 4.4, 0.0, 1.0),
                        masks=({"fit": (5000., 5001.), "synth": (5000., 5001.)},))


def _stub_result(resolved):
    from tumult.synthesis.request import SynthesisResult
    xx = np.linspace(5000., 5001., 5)
    return SynthesisResult((Spectrum(xx, np.ones_like(xx)),), (), resolved.request, {})


class _Ticks:
    """Collects the document's next-tick callbacks so a test runs them deliberately."""

    def __init__(self, document, monkeypatch):
        import threading
        self.queued, self.arrived, self.total = [], threading.Event(), 0

        def add(callback):
            self.queued.append(callback)
            self.total += 1
            self.arrived.set()
        monkeypatch.setattr(document, "add_next_tick_callback", add)

    def wait(self, total, timeout=10):
        """Whether ``total`` callbacks have been scheduled so far, within ``timeout``."""
        import time
        deadline = time.monotonic() + timeout
        while self.total < total and time.monotonic() < deadline:
            self.arrived.wait(0.2)
            self.arrived.clear()
        return self.total >= total

    def run(self):
        while self.queued:
            self.queued.pop(0)()

    def settle(self, predicate, timeout=10):
        """Run callbacks as they arrive until ``predicate()`` holds; whether it did."""
        import time
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            self.run()
            if predicate():
                return True
            self.arrived.wait(0.2)
            self.arrived.clear()
        self.run()
        return predicate()


def _send_action(document, value):
    """What the page does when a sidebar button is clicked: write the action channel."""
    document.select_one({"name": "object-actions"}).value = value


def test_four_overlapping_workflow_submissions_settle_on_the_document(tmp_path, monkeypatch):
    """Each durable synthesis is submitted once to the WORKFLOWS pool; four at once all
    run concurrently (none waits on another queued job) and deliver to the document."""
    import threading
    from tumult.synthesis import workflows

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    barrier = threading.Barrier(4, timeout=10)

    def execute(resolved, cancel_token=None):
        barrier.wait()                      # deadlocks unless all four workers are live
        return _stub_result(resolved)

    monkeypatch.setattr(workflows, "execute_request", execute)
    delivered = []
    handles = [ctx.submit_workflow(f"synth-{i}", _stub_request(tmp_path), delivered.append)
               for i in range(4)]
    active = ctx.jobs.active
    assert len(active) == 4 and [job.handle for job in active] == handles
    assert [job.workflow_id for job in active] == [h.id for h in handles]
    assert len({h.id for h in handles}) == 4
    assert ticks.wait(4), "every workflow must settle and schedule its delivery"
    assert all(h.done and h.status == "finished" for h in handles)
    assert not delivered                    # nothing reaches the workspace off-thread
    ticks.run()
    assert len(delivered) == 4 and all(len(r.spectra) == 1 for r in delivered)
    assert not len(ctx.jobs)


def test_queued_workflow_cancelled_from_the_processes_panel_settles(tmp_path, monkeypatch):
    """Cancelling a still-queued workflow through its Processes row cancels its future,
    records it as cancelled, and clears the row and the live handle without a worker."""
    import threading
    from concurrent.futures import ThreadPoolExecutor
    from tumult.synthesis import submission, workflows

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    monkeypatch.setattr(workflows, "execute_request",
                        lambda *a, **k: pytest.fail("a cancelled workflow ran"))
    gate = threading.Event()
    pool = ThreadPoolExecutor(1)
    monkeypatch.setattr(submission, "WORKFLOWS", pool)
    blocker = pool.submit(gate.wait, 10)
    try:
        delivered = []
        handle = ctx.submit_workflow("queued", _stub_request(tmp_path), delivered.append)
        job, = ctx.jobs.active
        assert job.name == "queued" and handle.status == "queued"
        assert job.display_state() == "queued"
        _send_action(document, f"cancel-pending:{job.pid}:1")
        assert handle.done and handle.status == "cancelled"
        assert handle.id not in submission._live
        assert ticks.settle(lambda: not len(ctx.jobs))
        assert job.state == "cancelled" and not delivered
    finally:
        gate.set()
        blocker.result(10)
        pool.shutdown(wait=True)


def test_closing_a_session_leaves_a_durable_workflow_running(tmp_path, monkeypatch):
    import threading
    from tumult.synthesis import workflows

    document, ctx = _build_document()
    _Ticks(document, monkeypatch)
    started, release, tokens = threading.Event(), threading.Event(), []

    def execute(resolved, cancel_token=None):
        tokens.append(cancel_token)
        started.set()
        assert release.wait(10)
        return _stub_result(resolved)

    monkeypatch.setattr(workflows, "execute_request", execute)
    handle = ctx.submit_workflow("detached", _stub_request(tmp_path), lambda r: None)
    assert started.wait(10)
    for callback in list(document.session_destroyed_callbacks):
        callback(None)
    assert not tokens[0].cancelled
    release.set()
    assert handle.result(timeout=10).spectra


# --- the correction tool's chained stages ----------------------------------------------
def _correction_tool(ctx, monkeypatch, tmp_path):
    """The document's correction tool with its modal inputs and result rendering stubbed."""
    import abundcorrtool
    from tumult.abundance import AbundanceResult

    tool = ctx.tools._tools["abundcorr"]
    tool._pending_atom = SimpleNamespace(symbol="Mg", path=str(tmp_path / "mg.atom"))
    monkeypatch.setattr(tool.setup, "context", lambda: _stub_request(tmp_path))
    monkeypatch.setattr(tool.grids, "reference_grid", lambda: {"folder": str(tmp_path)})
    for name in ("add_spectrum", "set_step", "select_result", "_on_corr_done"):
        monkeypatch.setattr(tool, name, lambda *a, **k: None)
    calls = {"reference": [], "corrections": []}

    def measure(request, atom, *, ref_grid, opts=None, cancel_token=None):
        calls["reference"].append(opts)
        result = AbundanceResult(element="Mg", symbol="Mg", summary={}, uid="corr",
                                 atom_path=atom, profile=opts.profile)
        result.ref_spec = None              # the stage attaches the reference spectrum
        return result

    def correct(request, working, atom, *, opts=None, progress=None, cancel_token=None):
        calls["corrections"].append(opts)
        return working

    monkeypatch.setattr(abundcorrtool, "measure_reference_lines", measure)
    monkeypatch.setattr(abundcorrtool, "run_nlte_corrections", correct)
    return tool, calls


def test_correction_options_are_captured_when_the_modal_submits(tmp_path, monkeypatch):
    """Both stages run on the options captured at submission; editing the reopened modal
    while the reference runs cannot change the admitted correction."""
    import threading
    import abundcorrtool

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    tool, calls = _correction_tool(ctx, monkeypatch, tmp_path)
    tool.form1.set_values(profile="gauss", ew_min=7.0, limit=0.75, line_grid=False,
                          linelist=True)
    started, release = threading.Event(), threading.Event()
    measure = abundcorrtool.measure_reference_lines

    def slow_measure(*args, **kwargs):
        started.set()
        assert release.wait(10)
        return measure(*args, **kwargs)

    monkeypatch.setattr(abundcorrtool, "measure_reference_lines", slow_measure)
    tool.run()
    assert started.wait(10)                 # the reference is running
    tool.form1.set_values(limit=2.5, profile="trapz", ew_min=1.0)
    release.set()
    assert ticks.settle(lambda: len(calls["corrections"]) == 1 and not len(ctx.jobs))
    for opts in calls["reference"] + calls["corrections"]:
        assert (opts.profile, opts.ew_min, opts.limit, opts.line_grid, opts.linelist) == (
            "gauss", 7.0, 0.75, False, True)


def test_cancelling_the_reference_stage_prevents_the_correction_stage(tmp_path, monkeypatch):
    import threading
    import abundcorrtool

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    tool, calls = _correction_tool(ctx, monkeypatch, tmp_path)
    started, release = threading.Event(), threading.Event()
    measure = abundcorrtool.measure_reference_lines

    def slow_measure(*args, **kwargs):
        started.set()
        assert release.wait(10)
        return measure(*args, **kwargs)     # a worker that never looks at its token

    monkeypatch.setattr(abundcorrtool, "run_nlte_corrections",
                        lambda *a, **k: pytest.fail("stage 2 admitted after cancellation"))
    monkeypatch.setattr(abundcorrtool, "measure_reference_lines", slow_measure)
    tool.run()
    assert started.wait(10)
    job, = ctx.jobs.active
    _send_action(document, f"cancel-pending:{job.pid}:1")
    assert job.state == "cancelling"        # the worker has not noticed yet
    release.set()
    assert ticks.settle(lambda: not len(ctx.jobs))
    assert job.state == "cancelled"
    assert len(calls["reference"]) == 1 and not calls["corrections"]
    assert not tool.results


# --- programmatic writes to the drag sources -------------------------------------------
def test_callback_suppression_restores_the_previous_flag():
    from filtertool import FilterTool
    from normalisetool import NormaliseTool
    from tool import suppress_callbacks

    class Failing:
        _suppress = False

        @suppress_callbacks
        def write(self):
            assert self._suppress
            raise ValueError("write failed")

        @suppress_callbacks
        def outer(self):
            self.inner()
            return self._suppress

        @suppress_callbacks
        def inner(self):
            return "inner"

    obj = Failing()
    with pytest.raises(ValueError):
        obj.write()
    assert obj._suppress is False           # a failed write cannot mute later drags
    assert obj.outer() is True and obj._suppress is False

    dragged = []
    filt = FilterTool.__new__(FilterTool)
    filt._suppress, filt.active = False, True
    filt.cross = ColumnDataSource(dict(x=[0.], y=[0.]))
    filt.vline, filt.hline = SimpleNamespace(location=0.), SimpleNamespace(location=0.)
    filt._debounce = lambda key, fn: dragged.append(key)
    filt.cross.on_change("data", filt._on_drag)
    filt._set_cross(1., 2.)
    assert not dragged and filt.cross.data["x"] == [1.]
    filt.cross.data = dict(x=[3.], y=[4.])  # a user drag still refits
    assert dragged == ["refit"]

    norm = NormaliseTool.__new__(NormaliseTool)
    norm._suppress, norm.active = False, True
    norm.knots = ColumnDataSource(dict(x=[], y=[]))
    norm._debounce = lambda key, fn: dragged.append(key)
    norm.knots.on_change("data", norm._on_knots)
    norm._set_knots([1., 2.], [1., 1.])
    assert dragged == ["refit"] and norm.knots.data["x"] == [1., 2.]
    norm.knots.data = dict(x=[1.], y=[1.])
    assert len(dragged) == 2 and norm._suppress is False


# --- job outcomes and diagnostics -------------------------------------------------------
def _footer(document):
    return document.select_one({"name": "footer-status"}).value


def _processes(document):
    return document.select_one({"name": "processes"}).value


def test_failed_loader_keeps_a_row_with_the_error_until_dismissed(monkeypatch, caplog):
    """A loader exception is logged with its traceback from the worker; the Processes
    row survives with the concise error, the footer points at it, and the row leaves
    on its dismiss button."""
    import logging
    from concurrent.futures import ThreadPoolExecutor

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    added = []

    def loader(token):
        raise ValueError("bad spectrum file")

    with caplog.at_level(logging.ERROR, logger="tumultstudio.jobs"):
        with ThreadPoolExecutor(1) as pool:
            future = ctx.submit("broken-load", loader, added.append)
            future.result(10)
        assert ticks.settle(lambda: not ctx.jobs.active)
    record, = [r for r in caplog.records if "broken-load" in r.getMessage()]
    assert record.exc_info and record.exc_info[0] is ValueError
    job, = ctx.jobs.retained
    assert job.state == "failed" and job.error == "ValueError: bad spectrum file"
    assert not added
    panel = _processes(document)
    assert "bad spectrum file" in panel and "dismiss-job" in panel
    assert "cancel-pending" not in panel and "spinner" not in panel
    assert "broken-load failed" in _footer(document) and "Processes" in _footer(document)

    _send_action(document, f"dismiss-job:{job.pid}:1")
    assert not len(ctx.jobs) and "dismiss-job" not in _processes(document)


def test_admission_failures_settle_as_failed_rows(monkeypatch, tmp_path):
    """A pool that refuses the job, or a workflow submission that raises, ends on the
    same terminal path as a failure in flight: no spinner, a retained row with the cause."""
    from concurrent.futures import ThreadPoolExecutor
    import app

    document, ctx = _build_document()
    _Ticks(document, monkeypatch)
    closed = ThreadPoolExecutor(1)
    closed.shutdown()
    monkeypatch.setattr(app, "_SYNTH_EXEC", closed)
    assert ctx.submit("refused", lambda token: None, lambda *_: None) is None
    job, = ctx.jobs
    assert job.state == "failed" and "RuntimeError" in job.error
    assert not ctx.jobs.active

    def refuse(request):
        raise OSError("workflow root is read-only")

    monkeypatch.setattr(app, "submit_synthesis", refuse)
    assert ctx.submit_workflow("unsubmittable", _stub_request(tmp_path), lambda r: None) is None
    _, job = ctx.jobs
    assert job.state == "failed" and job.error == "OSError: workflow root is read-only"
    assert "workflow root is read-only" in _processes(document)
    assert "see Processes" in _footer(document)


def test_cancellation_request_stays_visible_until_the_worker_settles(monkeypatch):
    """A worker that never looks at its token keeps its row, labelled as cancelling and
    without a cancel button, until it returns; only then is the job cancelled."""
    import threading

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    started, release = threading.Event(), threading.Event()
    delivered = []

    def native(token):
        started.set()
        assert release.wait(10)             # an uninterruptible call
        return "late value"

    ctx.submit("native", native, delivered.append)
    assert started.wait(10)
    assert ticks.settle(lambda: ctx.jobs.active[0].state == "running")
    job, = ctx.jobs.active
    _send_action(document, f"cancel-pending:{job.pid}:1")
    assert job.state == "cancelling" and job.token.cancelled
    panel = _processes(document)
    assert "cancelling" in panel and "spinner" in panel and "cancel-pending" not in panel
    assert "Cancelling native" in _footer(document)
    _send_action(document, f"cancel-pending:{job.pid}:2")    # a second click is a no-op
    assert job.state == "cancelling"
    release.set()
    assert ticks.settle(lambda: not len(ctx.jobs))
    assert job.state == "cancelled" and not delivered


def test_failure_after_document_destruction_is_still_logged(monkeypatch, caplog):
    """A detached job failing after its session closed logs the traceback and records
    the terminal state without touching the closed document."""
    import logging
    import threading

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    started, release = threading.Event(), threading.Event()

    def detached(token):
        started.set()
        assert release.wait(10)
        raise RuntimeError("dispatch exited with status 3: /runs/job.log")

    with caplog.at_level(logging.ERROR, logger="tumultstudio.jobs"):
        future = ctx.submit("survivor", detached, lambda *_: None, detached=True)
        assert started.wait(10)
        for callback in list(document.session_destroyed_callbacks):
            callback(None)
        assert ctx.closed
        scheduled = ticks.total
        release.set()
        future.result(10)                   # work() reports; it does not raise
    job, = ctx.jobs
    assert job.state == "failed" and "status 3" in job.error
    assert ticks.total == scheduled         # nothing scheduled on the closed document
    record, = [r for r in caplog.records if "survivor" in r.getMessage()]
    assert record.exc_info and record.exc_info[0] is RuntimeError


def test_errors_reach_only_their_own_document(monkeypatch):
    """Two documents on one server: a failure in one leaves the other's footer alone."""
    from concurrent.futures import ThreadPoolExecutor

    doc_a, ctx_a = _build_document()
    ticks_a = _Ticks(doc_a, monkeypatch)
    doc_b, ctx_b = _build_document()
    ticks_b = _Ticks(doc_b, monkeypatch)

    def loader(token):
        raise ValueError("only document A asked for this")

    with ThreadPoolExecutor(1):
        future = ctx_a.submit("private-failure", loader, lambda *_: None)
        future.result(10)
    assert ticks_a.settle(lambda: not ctx_a.jobs.active)
    assert "private-failure failed" in _footer(doc_a)
    assert ticks_b.total == 0 and not len(ctx_b.jobs)
    assert "private-failure" not in _footer(doc_b)
    assert not [h for h in __import__("logging").getLogger().handlers
                if type(h).__name__ == "FooterErrorHandler"]


def test_late_completion_cannot_mutate_a_superseded_workspace(monkeypatch):
    """Work admitted before a workspace boundary settles its row but delivers nothing;
    result-job progress published under the old generation is dropped as well."""
    import threading
    from jobs import ResultJobs

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    release, delivered = threading.Event(), []

    def slow(token):
        assert release.wait(10)
        return "stale"

    ctx.submit("stale-load", slow, delivered.append)
    old = ctx.generation
    ctx.advance_generation()                # what a session load does at its commit
    assert ctx.jobs.active[0].generation is old
    release.set()
    assert ticks.settle(lambda: not len(ctx.jobs))
    assert not delivered

    queued, callbacks, rendered, completed = [], [], [], []
    result = SimpleNamespace(uid="fit", lines=[SimpleNamespace(a=0.)], termination="measured")
    stub = SimpleNamespace(submit=lambda name, run, done, detached=False: queued.append((run, done)),
                           schedule=callbacks.append, notify=lambda *a: None, closed=False,
                           generation=object())
    jobs = ResultJobs(stub, lambda uid: result, lambda r: rendered.append(r.lines[0].a))

    def work(working, token, progress):
        working.lines[0].a = 1.
        progress(working)
        return working

    jobs.submit("fit", work, completed.append, result=result)
    run, done = queued.pop()
    outcome = run(CancelToken())
    stub.generation = object()              # the workspace was replaced meanwhile
    callbacks.pop()()
    done(outcome)
    assert result.lines[0].a == 0. and not rendered and not completed


# --- save slots: acknowledgements follow what happened on disk ----------------------
def _session_action(document, value):
    document.select_one({"name": "session-actions"}).value = value


def _session_list(document):
    return document.select_one({"name": "session-render"}).value


def _with_spectrum(ctx):
    x = np.arange(5000., 5010., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="obs",
                          source=("file", "/nowhere/obs.fits"))


def test_save_is_acknowledged_only_after_the_config_is_replaced(tmp_path, monkeypatch):
    """A failed replacement reports the destination and cause, leaves the previous config
    and the live workspace untouched, and the retry after the cause is fixed succeeds."""
    import store

    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    _with_spectrum(ctx)
    _session_action(document, "save first 1")
    assert "Saved “first”" in _footer(document)
    assert 'data-sid="1"' in _session_list(document)
    before = store.CONFIG_PATH.read_text()
    replace = store.os.replace

    def fail(src, dst):
        raise OSError(28, "No space left on device")
    monkeypatch.setattr(store.os, "replace", fail)
    _session_action(document, "save second try 2")
    footer = _footer(document)
    assert "Save “second try” failed" in footer and "No space left on device" in footer
    assert str(store.CONFIG_PATH) in footer
    assert store.CONFIG_PATH.read_text() == before
    assert 'data-sid="2"' not in _session_list(document)
    assert [o.name for o in ctx.canvas.objects] == ["obs"]
    assert not list(tmp_path.glob("*.tmp"))

    monkeypatch.setattr(store.os, "replace", replace)  # the disk has room again
    _session_action(document, "save second try 3")
    assert "Saved “second try”" in _footer(document)
    assert 'data-sid="2"' in _session_list(document)
    assert [s["name"] for s in store.list_sessions()] == ["second try", "first"]


def test_delete_failure_is_visible_and_keeps_the_slot(tmp_path, monkeypatch):
    import store

    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    _session_action(document, "save keep 1")
    before = store.CONFIG_PATH.read_text()

    def fail(src, dst):
        raise PermissionError(13, "Permission denied")
    monkeypatch.setattr(store.os, "replace", fail)
    _session_action(document, "delete 1 2")
    footer = _footer(document)
    assert "Delete failed" in footer and "Permission denied" in footer
    assert store.CONFIG_PATH.read_text() == before
    assert 'data-sid="1"' in _session_list(document)
    assert store.get_session(1) is not None


def test_broken_config_is_reported_and_never_replaced_from_the_gui(tmp_path, monkeypatch):
    """A config the store rejects (here an unsupported version) blocks saving, loading,
    deleting and refreshing with the same named cause; the file is left as it was."""
    import store

    path = tmp_path / "cfg.json"
    monkeypatch.setattr(store, "CONFIG_PATH", path)
    document, ctx = _build_document()
    _with_spectrum(ctx)
    path.write_text('{"version": 7, "preferences": {}, "recent_files": [], "sessions": []}')
    for action in ("save x 1", "load 1 2", "delete 1 3", "refresh 4"):
        _session_action(document, action)
        footer = _footer(document)
        assert "version 7 is not supported" in footer and str(path) in footer, action
    assert '"version": 7' in path.read_text()
    assert [o.name for o in ctx.canvas.objects] == ["obs"]
