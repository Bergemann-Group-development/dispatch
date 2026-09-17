"""Background computations must not apply stale UI state or wait on HOP locks."""
import html
import itertools
import os
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from bokeh.models import ColumnDataSource, Range1d, TapTool

from tumult import Spectrum
from tumult.synthesis import AtmosphereGrid, SynthConfig, SynthRequest
from tumult.synthesis.cancel import CancelToken


@pytest.fixture(autouse=True)
def _temporary_config(tmp_path, monkeypatch):
    """Every document built here reads and writes a config of its own; the project's
    ``.tumultstudio.json`` is neither consulted nor polluted."""
    import store
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "tumultstudio.json")


@pytest.fixture(autouse=True)
def _daemon_in_this_process(workflow_daemon, monkeypatch):
    """Workflows run on a daemon thread of the test process, where the stubbed
    calculations live. Rows are drawn as soon as a job is admitted; the delay that hides
    jobs settling at once has its own test."""
    import processes
    monkeypatch.setattr(processes, "SHOW_AFTER_S", 0.0)
    return workflow_daemon


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
    grids.reference_grid = lambda: None
    params = {"teff": 5750.0, "logg": 4.5, "feh": 0.0, "vmic": 1.0}
    monkeypatch.setattr(ionbalance, "synth_config", lambda prefs: SynthConfig())
    context = setup.context(params)

    assert setup.form.value("atmos_grid_id") == "stagger"
    assert [value for value, _ in setup.form.fields["atmos_grid_id"].options] == [
        "marcs", "stagger"]
    # an interpolated grid, LTE and scaled-solar, starting at the entry point
    assert isinstance(context.atmosphere, AtmosphereGrid)
    assert context.atmosphere.folder == Path("/stagger")
    assert context.departures is None and context.abundances is None
    assert dict(context.star) == params
    # background lines are off unless switched on: the ladders hold the element's lines only
    assert setup.background_lines is False
    setup.form.fields["background_lines"].value = True
    assert setup.background_lines is True


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
    # the job's windows keep each mask's fit and synth bounds apart; a fit synthesises the
    # span of the synth windows and compares inside the fit windows
    assert list(modal._job_windows()) == [
        ((5004.5, 5005.5), (5004.0, 5006.0)), ((5001.0, 5002.0), (5000.5, 5002.5))]


# --- durable synthesis workflows on the document ---------------------------------------
def _stub_request(root, windows=((5000., 5001.),), **preferences):
    """A synthesis request whose operation the tests stub out."""
    return SynthRequest(config=SynthConfig(root=root, preferences=preferences),
                        atmosphere=AtmosphereGrid("grid"), star=(5750.0, 4.4, 0.0, 1.0),
                        windows=windows)


def _stub_result(request):
    from tumult.synthesis.request import SynthesisResult
    xx = np.linspace(5000., 5001., 5)
    return SynthesisResult((Spectrum(xx, np.ones_like(xx)),), request.windows, request, {}, {})


def _submit_synthesis(ctx, name, request, on_done, meta=None):
    """Submit a synthesis through the document; its workflow handle, or ``None``."""
    job = ctx.submit_operation(name, "synthesis.synthesize", request, on_done, meta=meta)
    return None if job is None else job.handle


def _stub_synthesize(monkeypatch, synthesize):
    from tumult.synthesis import operations
    monkeypatch.setattr(operations, "synthesize", synthesize)


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


def test_overlapping_workflow_submissions_settle_on_the_document(tmp_path, monkeypatch):
    """Four durable syntheses submitted at once each run in the daemon, one after the
    other, and deliver their result on the document thread only, once each."""

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    order = []

    def execute(request, background_lines=True, cancel_token=None):
        order.append(request)
        return _stub_result(request)

    _stub_synthesize(monkeypatch, execute)
    delivered = []
    handles = [_submit_synthesis(ctx, f"synth-{i}", _stub_request(tmp_path), delivered.append)
               for i in range(4)]
    active = ctx.jobs.active
    assert len(active) == 4 and [job.handle for job in active] == handles
    assert [job.workflow_id for job in active] == [h.id for h in handles]
    assert len({h.id for h in handles}) == 4
    for handle in handles:
        handle.result(timeout=10)
    assert all(h.status == "finished" for h in handles) and len(order) == 4
    assert not delivered                    # nothing reaches the workspace off-thread
    assert ticks.settle(lambda: len(delivered) == 4)
    assert all(len(r.spectra) == 1 for r in delivered)
    assert ticks.settle(lambda: not len(ctx.jobs))


def test_queued_workflow_cancelled_from_the_processes_panel_settles(tmp_path, monkeypatch):
    """A workflow queued behind a running one is cancelled from its Processes row: the
    daemon settles its record as cancelled without running it, and the row clears."""
    import threading

    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    started, release = threading.Event(), threading.Event()

    def execute(request, background_lines=True, cancel_token=None):
        assert request.config.preferences.get("first"), "a cancelled workflow ran"
        started.set()
        assert release.wait(10)
        return _stub_result(request)

    _stub_synthesize(monkeypatch, execute)
    first = _stub_request(tmp_path, first=True)
    _submit_synthesis(ctx, "running", first, lambda r: None)
    assert started.wait(10)
    delivered = []
    handle = _submit_synthesis(ctx, "queued", _stub_request(tmp_path), delivered.append)
    job = ctx.jobs.find_workflow(handle.id)
    assert job.name == "queued" and handle.status == "queued"
    assert job.display_state() == "queued"
    assert "Synthesis · queued" in _processes(document)
    _send_action(document, f"cancel-pending:{job.pid}:1")
    assert handle.done and handle.status == "cancelled"
    assert ticks.settle(lambda: job.state == "cancelled")
    assert not delivered
    release.set()
    assert ticks.settle(lambda: not ctx.jobs.active)


def test_closing_a_session_leaves_a_durable_workflow_running(tmp_path, monkeypatch):
    import threading

    document, ctx = _build_document()
    _Ticks(document, monkeypatch)
    started, release, tokens = threading.Event(), threading.Event(), []

    def execute(request, background_lines=True, cancel_token=None):
        tokens.append(cancel_token)
        started.set()
        assert release.wait(10)
        return _stub_result(request)

    _stub_synthesize(monkeypatch, execute)
    handle = _submit_synthesis(ctx, "detached", _stub_request(tmp_path), lambda r: None)
    assert started.wait(10)
    for callback in list(document.session_destroyed_callbacks):
        callback(None)
    assert not tokens[0].cancelled
    release.set()
    assert handle.result(timeout=10).spectra


# --- the correction tool's chained stages ----------------------------------------------
def _stub_corrections_inputs(monkeypatch, tool, tmp_path):
    """The correction tool's request without windows, a reference grid, and neither a
    preflight nor a quadrature preparation, which need files on disk."""
    from tumult.departures import DepartureRequest
    monkeypatch.setattr(tool.setup, "context", lambda: _stub_request(tmp_path, windows=()))
    monkeypatch.setattr(tool.grids, "reference_grid", lambda: {"folder": str(tmp_path)})
    monkeypatch.setattr(SynthRequest, "check", lambda self: self)
    monkeypatch.setattr(DepartureRequest, "prepare", lambda self, request: (None, None))


def _correction_tool(ctx, monkeypatch, tmp_path):
    """The document's correction tool with its modal inputs and result rendering stubbed.
    The stages are patched where the task runs them (``tumult.abundance``); the record
    runs in this process under the test preferences."""
    from tumult import abundance as abundcorrtool
    from tumult.abundance import AbundanceResult

    tool = ctx.tools._tools["abundcorr"]
    tool._pending_atom = SimpleNamespace(symbol="Mg", path=str(tmp_path / "mg.atom"))
    _stub_corrections_inputs(monkeypatch, tool, tmp_path)
    for name in ("add_spectrum", "set_step", "select_result", "_on_corr_done"):
        monkeypatch.setattr(tool, name, lambda *a, **k: None)
    calls = {"reference": [], "corrections": []}

    def measure(request, atom, *, opts=None, quadrature=None, cancel_token=None):
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
    from tumult import abundance as abundcorrtool

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
    from tumult import abundance as abundcorrtool

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


def _usage(document):
    return document.select_one({"name": "processes-usage"}).value


def _details(document):
    return document.select_one({"name": "job-details"}).value


def _recovery(document):
    return document.select_one({"name": "recovery-list"}).value


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

    def refuse(operation, request, **kwargs):
        raise OSError("workflow root is read-only")

    monkeypatch.setattr(app, "submit_operation_record", refuse)
    assert _submit_synthesis(ctx, "unsubmittable", _stub_request(tmp_path), lambda r: None) is None
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


def _session_load(document, sid, nonce):
    """Load slot ``sid`` the way the page does: a workspace with unsaved changes gets the
    Save / Discard / Cancel prompt first, and these tests always discard."""
    _session_action(document, f"load {sid} {nonce}")
    _session_action(document, f"load-discard {nonce}")


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
        if action.startswith("load"):
            _session_load(document, *action.split()[1:])
        else:
            _session_action(document, action)
        footer = _footer(document)
        assert "version 7 is not supported" in footer and str(path) in footer, action
    assert '"version": 7' in path.read_text()
    assert [o.name for o in ctx.canvas.objects] == ["obs"]


def test_a_rejected_config_does_not_block_opening_a_file(tmp_path, monkeypatch):
    """The recent-files entry is a convenience: when the store refuses the config, the
    footer says so at warn level, the config is left alone and the pick still loads."""
    import store

    path = tmp_path / "cfg.json"
    monkeypatch.setattr(store, "CONFIG_PATH", path)
    document, ctx = _build_document()
    ticks = _Ticks(document, monkeypatch)
    path.write_text('{"version": 7, "preferences": {}, "recent_files": [], "sessions": []}')
    spectrum = _spectrum_file(tmp_path)
    picker = document.select_one({"name": "file_picker"})
    picker.children[2].value = f"{spectrum}\x001"
    assert ticks.settle(lambda: [o.name for o in ctx.canvas.objects] == ["obs.txt"])
    footer = _footer(document)
    assert "Recent files not updated" in footer and "version 7 is not supported" in footer
    assert "status-warn" in footer
    assert '"recent_files": []' in path.read_text()
    assert not ctx.jobs.retained


# --- save slots: the whole object workspace, restored transactionally -----------------
ATOM_FIXTURE = Path(__file__).resolve().parents[1] / "input" / "atoms" / "atom.ba06"


def _spectrum_file(tmp_path, name="obs.txt", n=400):
    path = tmp_path / name
    path.write_text("".join(f"{5000 + 0.01 * i:.3f} {1.0 - 0.3 * (i % 50 == 25):.3f}\n"
                            for i in range(n)))
    return path


def _mask_file(tmp_path, name="fe-lines.txt"):
    path = tmp_path / name
    path.write_text("; element: Fe\n5001.0 5000.5 5001.5 Fe I\n5002.0 5001.8 5002.2 Fe II\n")
    return path


def _diag_line(lo=5000.5, hi=5001.5):
    from tumult.abundance import LineResult
    from tumult.linemask import LineMask
    mask = LineMask(lam0=0.5 * (lo + hi))
    mask.set_fit(lo, hi, 0.5)
    return LineResult(lam0=0.5 * (lo + hi), ion="Fe I", e_low=1.0, loggf=-1.0, ew=50.0,
                      mask=mask)


def _populate(ctx, tmp_path):
    """Fill a document with one object of every kind, edited: a file spectrum, a synth
    spectrum (knots, convolution, hidden), an imported mask set edited away from its
    file, a model atom, and one result per tool; the correction result is focused on its
    Lines step with its first line jumped to. Returns the identities to check later."""
    from tumult.abundance import AbundanceResult
    from tumult.atom import AtomClass
    from tumult.measurement import BalanceResult
    from tumult.spectrum_io import load_file

    canvas, tools = ctx.canvas, ctx.tools
    obs_path = _spectrum_file(tmp_path)
    obs = canvas.add_object(load_file(obs_path), name="obs", source=("file", str(obs_path)))
    canvas.set_style(obs, color="#ff8800", width=2.0)
    x = np.linspace(5000., 5004., 401)
    synth = canvas.add_object(Spectrum(x, 2.0 * np.ones_like(x)), name="model",
                              source=("synth", None), meta={"teff": 5400})
    synth.norm_knots = (np.array([5000.0, 5004.0]), np.array([2.0, 2.0]))
    canvas.normalise_object(synth)
    synth.conv = {"vshift": 1.5}
    canvas._refresh_object(synth)
    canvas.toggle_object(synth.uid)                  # hidden

    mset = ctx.masks.load_file(_mask_file(tmp_path))
    ctx.masks.set_fit(0, 5000.4, 5001.6)             # edited: the set reads as modified
    assert mset.dirty
    atomtool = tools._tools["atom"]
    atom = atomtool.add(str(ATOM_FIXTURE), AtomClass(str(ATOM_FIXTURE)), select=False)

    abund, corr, ionbal = (tools._tools[k] for k in ("abundance", "abundcorr", "ionbal"))
    res_a = AbundanceResult(element="Fe", symbol="Fe", summary={}, uid="res-a", name="Fe A",
                            ref_uid=obs.uid, ref_name="obs", lines=[_diag_line()])
    res_b = BalanceResult(element="Fe", uid="res-b", name="Fe balance", ref_uid=obs.uid,
                          ref_name="obs", maskset_uid=str(mset.uid), maskset_name=mset.name,
                          lines=[_diag_line()])
    res_c = AbundanceResult(element="Ba", symbol="Ba", summary={}, uid="res-c", name="Ba corr",
                            atom_path=str(ATOM_FIXTURE), lines=[_diag_line(), _diag_line(5002., 5003.)])
    abund.results.append(res_a)
    abund.selected_uid = res_a.uid
    ionbal.results.append(res_b)
    corr.results.append(res_c)
    corr.select_result(res_c.uid)
    corr.set_step(1)
    corr._jump_line("1")
    canvas.main.x_range.start, canvas.main.x_range.end = 5001.0, 5002.5
    return SimpleNamespace(obs=obs.uid, synth=synth.uid, mset=mset.uid, atom=atom.uid,
                           obs_path=str(obs_path))


def test_object_panel_close_actions_leave_each_view(tmp_path, monkeypatch):
    """The × of every object view (mask set, atom, the three result views) reaches the
    ``object-actions`` channel as the tool's name and stands the tool down, leaving the
    sidebar without an active tool and the mode class off."""
    document, ctx = _build_document()
    ids = _populate(ctx, tmp_path)
    tools = ctx.tools
    tools._tools["atom"].select(ids.atom)
    tools._tools["ionbal"].select_result(tools._tools["ionbal"].results[0].uid)
    for name in ("linemask", "atom", "abundance", "abundcorr", "ionbal"):
        tool = tools._tools[name]
        tools.activate(name)
        assert tools.active_name() == name, name
        assert all(b.value.startswith("on:") for b in tool._mode_bridges), name
        _send_action(document, f"{name}:-:{name}")
        assert tools.active_name() is None and not tool.active, name
        assert all(b.value.startswith("off:") for b in tool._mode_bridges), name


def _saved_slot(tmp_path, monkeypatch, populate=_populate):
    """Save a populated document into a fresh config; return the ids and the slot."""
    import store
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    ids = populate(ctx, tmp_path)
    _session_action(document, "save mixed 1")
    assert "Saved “mixed”" in _footer(document)
    return ids, store.get_session(1)


def _state(ctx):
    """What a restored document must reproduce, keyed for comparison."""
    canvas, tools = ctx.canvas, ctx.tools
    abund, corr, ionbal = (tools._tools[k] for k in ("abundance", "abundcorr", "ionbal"))
    return {
        "objects": [(o.uid, o.name, o.visible, o.line.glyph.line_color, o.line.glyph.line_width,
                     o.conv, o.already_normalised,
                     None if o.norm_knots is None else np.round(o.norm_knots[0], 6).tolist())
                    for o in canvas.objects],
        "selected": canvas.selected_uid,
        "masks": [(s.uid, s.name, s.element, s.dirty, [(m.lmin, m.lmax) for m in s.masks])
                  for s in ctx.masks.listed],
        "mask_current": (ctx.masks.current_uid, ctx.masks.active),
        "atoms": [(a.uid, a.path) for a in tools._tools["atom"].atoms],
        "results": {name: ([(r.uid, r.ref_uid) for r in w.results], w.selected_uid, w.step)
                    for name, w in (("abundance", abund), ("abundcorr", corr),
                                    ("ionbal", ionbal))},
        "balance_mask": [r.maskset_uid for r in ionbal.results],
        "focus": (ctx.focus_kind, tools.active_name()),
        "jumped": corr._jumped,
        "main_x": (round(canvas.main.x_range.start, 3), round(canvas.main.x_range.end, 3)),
        "units": (canvas.xlabel, canvas.ylabel),
    }


@pytest.mark.skipif(not ATOM_FIXTURE.is_file(), reason="input/atoms/atom.ba06 is absent")
def test_mixed_workspace_round_trips_into_a_fresh_document(tmp_path, monkeypatch):
    """Every collection, its edits, the selections, the focused view and the explicit
    range come back in a fresh document; identities are the saved ones."""
    ids, slot = _saved_slot(tmp_path, monkeypatch)
    assert slot["schema_version"] == 2 and slot["focus"] == {"kind": "corr", "uid": "res-c"}
    assert slot["tools"]["jumped"]["abundcorr"] == 1

    document, ctx = _build_document()
    expected_ids = ids
    before = ctx.generation
    _session_load(document, 1, 2)
    footer = _footer(document)
    assert "Loaded “mixed”: 2 spectra, 1 line mask, 1 atom, 3 results" in footer
    assert "skipped" not in footer and ctx.generation is not before
    state = _state(ctx)
    assert [o[0] for o in state["objects"]] == [expected_ids.obs, expected_ids.synth]
    assert state["objects"][0][2:5] == (True, "#ff8800", 2.0)
    assert state["objects"][1][2] is False and state["objects"][1][5] == {"vshift": 1.5}
    assert state["objects"][1][7] == [5000.0, 5004.0]
    assert state["masks"] == [(expected_ids.mset, "fe-lines", "Fe", True,
                               [(5000.4, 5001.6), (5001.8, 5002.2)])]
    assert state["mask_current"] == (expected_ids.mset, 0)
    assert state["atoms"] == [(expected_ids.atom, str(ATOM_FIXTURE))]
    assert state["results"] == {"abundance": ([("res-a", expected_ids.obs)], "res-a", 1),
                                "abundcorr": ([("res-c", None)], "res-c", 1),
                                "ionbal": ([("res-b", expected_ids.obs)], None, 1)}
    assert state["balance_mask"] == [str(expected_ids.mset)]
    assert state["focus"] == ("corr", "abundcorr") and state["jumped"] == 1
    assert state["main_x"] == (5001.0, 5002.5)
    # the correction is listed under its atom, and the reference spectrum relinked
    atom_objects = ctx.tools._tools["atom"].objects_html()
    assert 'data-atom-results="1"' in atom_objects and "Corrections" not in atom_objects
    assert ctx.tools._tools["ionbal"]._source_set(ctx.tools._tools["ionbal"].results[0]) \
        is ctx.masks.mset


@pytest.mark.skipif(not ATOM_FIXTURE.is_file(), reason="input/atoms/atom.ba06 is absent")
def test_loading_replaces_an_unrelated_workspace_and_an_empty_slot_clears_it(tmp_path,
                                                                             monkeypatch):
    import store
    ids, slot = _saved_slot(tmp_path, monkeypatch)
    document, ctx = _build_document()
    x = np.arange(5000., 5010., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="other",
                          source=("file", "/nowhere/other.fits"))
    ctx.masks.add_set("scratch", "Ni")
    ctx.masks.add(5003., 5004.)
    tools = ctx.tools
    tools._tools["abundance"].results.append(
        __import__("tumult.abundance", fromlist=["AbundanceResult"]).AbundanceResult(
            element="Ni", symbol="Ni", summary={}, uid="other-res", name="Ni"))
    tools.toggle("linemask")
    _session_action(document, "save other 2")
    assert tools.active_name() == "linemask" and ctx.focus_kind == "linemask"

    _session_load(document, 1, 3)
    state = _state(ctx)
    assert [o[1] for o in state["objects"]] == ["obs", "model"]
    assert [m[1] for m in state["masks"]] == ["fe-lines"]
    assert state["results"]["abundance"][0] == [("res-a", ids.obs)]
    assert state["focus"] == ("corr", "abundcorr")

    # the empty slot: every collection replaced by nothing, the plain canvas back
    store.save_session("empty", {"schema_version": 2})
    _session_load(document, 3, 4)
    assert "Loaded “empty”: 0 spectra, 0 line masks, 0 atoms, 0 results" in _footer(document)
    state = _state(ctx)
    assert state["objects"] == [] and state["masks"] == [] and state["atoms"] == []
    assert all(r[0] == [] and r[1] is None for r in state["results"].values())
    assert state["focus"] == ("spectrum", None) and state["selected"] is None
    # and the other way round: loading the "other" slot into the emptied document
    _session_load(document, 2, 5)
    assert "skipped spectrum other (/nowhere/other.fits)" in _footer(document)
    state = _state(ctx)
    assert state["objects"] == [] and [m[1:3] for m in state["masks"]] == [("scratch", "Ni")]
    assert state["results"]["abundance"][0] == [("other-res", None)]
    assert state["focus"] == ("linemask", "linemask")


@pytest.mark.skipif(not ATOM_FIXTURE.is_file(), reason="input/atoms/atom.ba06 is absent")
def test_missing_files_skip_their_objects_and_clear_what_pointed_at_them(tmp_path,
                                                                         monkeypatch):
    """The partial-load policy: a spectrum or atom file that is gone is skipped and named
    with its path and reason; independent objects load; references and the focus that
    named a skipped object are cleared, never redirected."""
    import store

    def populate(ctx, tmp_path):
        ids = _populate(ctx, tmp_path)
        ctx.tools._tools["range"]  # noqa: the focus goes to the file spectrum
        ctx.canvas.select_object(ids.obs)
        ctx.tools.deactivate_active()
        ctx.tools._tools["abundcorr"].blur()
        ctx.focus_kind = "spectrum"
        return ids
    ids, slot = _saved_slot(tmp_path, monkeypatch, populate)
    assert slot["focus"] == {"kind": "spectrum", "uid": ids.obs}
    Path(ids.obs_path).unlink()
    slot["atoms"][0]["path"] = str(tmp_path / "no-such-atom")
    store.save_session("broken links", slot)

    document, ctx = _build_document()
    _session_load(document, 2, 2)
    footer = html.unescape(_footer(document))
    assert "1 spectrum, 1 line mask, 0 atoms, 3 results" in footer
    assert f"skipped spectrum obs ({ids.obs_path}): FileNotFoundError" in footer
    assert f"atom no-such-atom ({tmp_path / 'no-such-atom'}): FileNotFoundError" in footer
    assert "3 references cleared" in footer
    state = _state(ctx)
    assert [o[1] for o in state["objects"]] == ["model"] and state["selected"] is None
    assert state["focus"] == ("spectrum", None)
    assert state["results"]["abundance"][0] == [("res-a", None)]
    assert state["results"]["ionbal"][0] == [("res-b", None)]
    assert state["balance_mask"] == [str(ids.mset)]           # the mask set did load
    assert "Corrections" in ctx.tools._tools["atom"].objects_html()   # orphaned, still listed

    # every spectrum missing: masks, atoms and results still load
    slot["spectra"][1]["kind"], slot["spectra"][1]["path"] = "file", "/nowhere/model.txt"
    slot["spectra"][1]["data"] = None
    slot["atoms"][0]["path"] = str(ATOM_FIXTURE)
    store.save_session("no spectra", slot)
    _session_load(document, 3, 3)
    state = _state(ctx)
    assert state["objects"] == [] and [m[1] for m in state["masks"]] == ["fe-lines"]
    assert [a[1] for a in state["atoms"]] == [str(ATOM_FIXTURE)]
    assert state["results"]["abundcorr"][0] == [("res-c", None)]
    assert 'data-atom-results="1"' in ctx.tools._tools["atom"].objects_html()


def test_broken_slots_are_refused_before_the_workspace_changes(tmp_path, monkeypatch):
    import store

    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    _with_spectrum(ctx)
    ctx.masks.add_set("kept", "Fe")
    before = _state(ctx)
    generation = ctx.generation
    store.save_session("old format", {"schema_version": 1, "objects": []})
    store.save_session("bad knots", {"schema_version": 2, "spectra": [
        {"uid": "s", "kind": "synth", "name": "m",
         "data": Spectrum(np.linspace(1., 2., 5), np.ones(5)).to_dict(),
         "norm_knots": {"x": [1.0, 2.0], "y": [1.0]}}]})
    store.save_session("bad result", {"schema_version": 2, "results": {"ionbal": {
        "results": [{"uid": "r", "schema_version": 1, "lines": ["not a line"]}]}}})
    store.save_session("dangling focus", {"schema_version": 2,
                                          "focus": {"kind": "result", "uid": "ghost"}})
    for sid, cause in ((1, "version 1 is not supported"), (2, "norm_knots"),
                       (3, "results.ionbal"), (4, "focus: result 'ghost'")):
        _session_load(document, sid, sid)
        footer = html.unescape(_footer(document))
        assert "Could not load" in footer and cause in footer, sid
        assert _state(ctx) == before and ctx.generation is generation
    _session_action(document, "refresh 5")
    assert "unsupported format (version 1)" in _session_list(document)


def test_failed_application_restores_the_previous_workspace(tmp_path, monkeypatch):
    """When applying a prepared workspace fails part-way, the previous one is put back
    (its objects, masks, results and focus); the boundary itself (generation, cancelled
    jobs) stands, since the old work's outcome had already been disowned."""
    import store

    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    ctx.canvas.add_object(Spectrum(np.linspace(5000., 5001., 50), np.ones(50)), name="one",
                          source=("synth", None))
    _session_action(document, "save lowpass 1")
    slot = store.get_session(1)
    slot["filter"] = {"active": True, "method": "lowpass", "res": 10000.0}
    store.save_session("with filter", slot)

    ctx.canvas.remove_object(ctx.canvas.objects[0].uid)
    kept = ctx.canvas.add_object(Spectrum(np.linspace(5002., 5003., 50), np.ones(50)),
                                 name="kept", source=("synth", None))
    ctx.masks.add_set("kept masks", "Ni")
    from tumult.abundance import AbundanceResult
    abund = ctx.tools._tools["abundance"]
    abund.results.append(AbundanceResult(element="Ni", symbol="Ni", summary={}, uid="kept-r",
                                         name="Ni", ref_uid=kept.uid))
    abund.select_result("kept-r")
    before = _state(ctx)
    filtertool = ctx.tools._tools["filter"]

    def explode(res):
        raise RuntimeError("no filter today")
    monkeypatch.setattr(filtertool, "apply_saved_lowpass", explode)
    generation = ctx.generation
    _session_load(document, 2, 2)
    footer = _footer(document)
    assert "Loading “with filter” failed: no filter today" in footer
    assert "the previous workspace was restored" in footer
    assert _state(ctx) == before and ctx.generation is not generation
    assert [o.uid for o in ctx.canvas.objects] == [kept.uid]


def test_a_worker_finishing_after_a_load_adds_nothing_to_the_new_workspace(tmp_path,
                                                                            monkeypatch):
    """Loading a slot is a workspace boundary: the document's pending jobs are cancelled,
    a late result is dropped, and the row settles as cancelled."""
    import threading
    import store
    from jobs import CANCELLED
    from tumult.spectrum_io import load_file

    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "cfg.json")
    document, ctx = _build_document()
    path = _spectrum_file(tmp_path)
    ctx.canvas.add_object(load_file(path), name="obs", source=("file", str(path)))
    _session_action(document, "save one 1")
    ticks = _Ticks(document, monkeypatch)
    release, delivered = threading.Event(), []

    def slow(token):
        assert release.wait(10)
        return Spectrum(np.linspace(5000., 5001., 20), np.ones(20))

    ctx.submit("late-load", slow,
               lambda spec: delivered.append(ctx.canvas.add_object(spec, name="late")))
    job = ctx.jobs.active[0]
    ctx.restore_workspace(store.get_session(1), "one")
    assert job.token.cancelled
    release.set()
    assert ticks.settle(lambda: job.state == CANCELLED)
    assert not delivered and [o.name for o in ctx.canvas.objects] == ["obs"]


# --- recovering detached workflows after a reconnect -----------------------------------
def _reconnecting_document(tmp_path, monkeypatch):
    """A fresh document whose configured synthesis root is ``tmp_path/synth``, where the
    stub requests of this file write their workflow records."""
    import store
    store.set_section("synthesis", {"dispatch_cwd": str(tmp_path)})
    document, ctx = _build_document()
    return document, ctx, _Ticks(document, monkeypatch)


def _blocking_execute(monkeypatch):
    """Stub execution that waits for ``release`` and honours its token; returns the
    events and the tokens it saw."""
    import threading
    started, release, tokens = threading.Event(), threading.Event(), []

    def execute(request, background_lines=True, cancel_token=None):
        tokens.append(cancel_token)
        started.set()
        assert release.wait(10)
        cancel_token.check()
        return _stub_result(request)

    _stub_synthesize(monkeypatch, execute)
    return started, release, tokens


def _destroy(document):
    for callback in list(document.session_destroyed_callbacks):
        callback(None)


def test_a_running_workflow_is_recovered_by_a_new_document_and_cancelled(tmp_path, monkeypatch):
    """The tab that submitted a workflow goes away; a reconnecting document finds the
    live record, attaches to it (a job of its own, cancellable through the handle) and
    its cancellation reaches the worker's token. The old document delivers nothing."""
    from jobs import CANCELLED

    started, release, tokens = _blocking_execute(monkeypatch)
    old_document, old_ctx, _ = _reconnecting_document(tmp_path, monkeypatch)
    old_delivered = []
    handle = _submit_synthesis(old_ctx, "first", _stub_request(tmp_path), old_delivered.append)
    assert started.wait(10)
    assert old_ctx.workflow_refs.to_list() == [
        {"id": handle.id, "folder": str(handle.folder), "name": "first", "imported": False}]
    _destroy(old_document)

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    panel = _processes(document)       # a row from the first paint: the daemon runs it
    assert "Synthesis · first" in panel and f'data-action="workflow-cancel" data-uid="{handle.id}"' in panel
    assert not ctx.jobs.active and not len(ctx.recovery)
    _send_action(document, f"job-details:{handle.id}:1")
    details = _details(document)
    assert "Synthesis · first" in details and 'data-action="workflow-attach"' in details
    _send_action(document, f"workflow-attach:{handle.id}:1")
    job, = ctx.jobs.active
    assert job.workflow_id == handle.id and job.generation is ctx.generation
    assert not len(ctx.recovery) and handle.id in ctx.workflow_refs
    assert "running" in _processes(document) and "cancel-pending" in _processes(document)

    _send_action(document, f"cancel-pending:{job.pid}:1")
    assert tokens[0].cancelled
    release.set()
    assert ticks.settle(lambda: job.state == CANCELLED)
    assert handle.status == "cancelled" and not old_delivered
    assert not ctx.canvas.objects


def test_a_running_correction_is_resumed_by_a_new_document_with_its_progress(tmp_path,
                                                                              monkeypatch):
    """The tab that started NLTE corrections goes away while the record runs; a new
    document attaches to it, the correction tool lists the record's stage-1 snapshot as
    a result (its reference lines), streams the later snapshots into it and receives
    the finished result through the tool's own delivery."""
    import threading
    from tumult import abundance
    from tumult.abundance import AbundanceResult

    store_root = tmp_path
    stage1, release, seen = threading.Event(), threading.Event(), []

    def measure(request, atom, *, opts=None, quadrature=None, cancel_token=None):
        return AbundanceResult(element="Mg", symbol="Mg", summary={}, uid="corr-1",
                               atom_path=atom, profile=opts.profile, stage="lines")

    def correct(request, working, atom, *, opts=None, progress=None, cancel_token=None):
        stage1.set()
        assert release.wait(10)
        working.stage = "corr"
        working.summary = {"a_mean": 0.1, "n": 0}
        return working

    monkeypatch.setattr(abundance, "measure_reference_lines", measure)
    monkeypatch.setattr(abundance, "run_nlte_corrections", correct)
    monkeypatch.setattr("tumult.synthesis.submission.PROGRESS_MIN_S", 0.0)

    old_document, old_ctx, _ = _reconnecting_document(store_root, monkeypatch)
    old_tool = old_ctx.tools._tools["abundcorr"]
    old_tool._pending_atom = SimpleNamespace(symbol="Mg", path=str(tmp_path / "mg.atom"))
    _stub_corrections_inputs(monkeypatch, old_tool, tmp_path)
    for name in ("add_spectrum", "set_step", "select_result"):
        monkeypatch.setattr(old_tool, name, lambda *a, **k: None)
    old_tool.run()
    job, = old_ctx.jobs.active
    handle = job.handle
    assert handle.kind.name == "abundance.nlte_corrections" and stage1.wait(10)
    _destroy(old_document)              # the old document never renders the outcome

    document, ctx, ticks = _reconnecting_document(store_root, monkeypatch)
    tool = ctx.tools._tools["abundcorr"]
    for name in ("add_spectrum", "set_step", "select_result"):
        monkeypatch.setattr(tool, name, lambda *a, **k: None)
    assert "NLTE corrections · Mg" in _processes(document)
    _send_action(document, f"workflow-attach:{handle.id}:1")
    job, = ctx.jobs.active
    assert job.workflow_id == handle.id and not tool.results
    # the stage-1 snapshot is on disk already: following the record reads it once, off
    # the document thread
    assert ticks.settle(lambda: bool(tool.results))
    res, = tool.results
    assert res.uid == "corr-1" and res.stage == "lines" and tool.selected_uid == "corr-1"

    release.set()
    assert ticks.settle(lambda: not ctx.jobs.active)
    assert tool.results[0] is res and res.stage == "corr"       # updated in place
    assert res.summary["a_mean"] == 0.1 and ctx.workflow_refs.imported(handle.id)
    assert handle.status == "finished" and not old_tool.results


def test_external_runs_are_listed_cancelled_and_reused_as_departures(tmp_path, monkeypatch,
                                                                    workflow_daemon):
    """A dispatch run started outside the GUI is found by the daemon through its own
    job.nml under the synthesis root and shown as a row: a live one can be ended through
    the daemon, a departure solve handed to the modal; a lost run and a run a workflow
    accounts for are not listed."""
    import json
    import socket
    from tumult.synthesis import daemon, jobrecord
    root = tmp_path / "synth"

    def native(folder, **fields):
        values = {"state": "running", "pid": os.getpid(), "host": socket.gethostname(),
                  "cwd": str(tmp_path), "input": f"{folder.name}.nml",
                  "started": "2026-09-13T11:15:20+02:00", "ended": "", "mpi_size": 1,
                  "omp_threads": 4, "slurm_job_id": "", "dispatch_hash": "abc", **fields}
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "job.nml").write_text("&job_record\n" + "".join(
            f" {k}={v!r},\n" if isinstance(v, str) else f" {k}={v},\n"
            for k, v in values.items()) + " /\n")

    native(root / "dep_fe607")                       # a running departure solve
    (root / "dep_fe607" / "save").mkdir()
    (root / "dep_fe607" / "save" / "pop_Fe_meta.txt").write_text("")
    native(root / "lte_5000", pid=2 ** 30)           # lost: its process is gone
    native(root / "owned")                           # named by a workflow record: not listed
    record = root / "_workflows" / "w1"
    record.mkdir(parents=True)
    (record / "workflow.json").write_text(json.dumps({
        "id": "w1", "state": "failed", "phase": "failed", "submitted": 1.0,
        "runs": [str(root / "owned")], "request": {}}))

    daemon.client().request("list", roots=[str(root)])
    workflow_daemon.scan_external()
    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    panel = _processes(document)
    run_id = jobrecord.run_id(root / "dep_fe607")
    assert "Dispatch run · dep_fe607" in panel and "4 threads" in panel
    assert "lte_5000" not in panel and "owned" not in panel
    assert f'data-action="run-cancel" data-uid="{run_id}"' in panel
    assert _usage(document).startswith("4/")        # outside the budget, still counted
    _send_action(document, f"job-details:{run_id}:1")
    assert 'data-action="run-departures"' in _details(document)

    _send_action(document, f"run-departures:{run_id}:1")
    setup = ctx.synth_modal.setup
    assert setup.mode_form.value("departures") == str(root / "dep_fe607")
    assert setup.is_nlte

    signalled = []
    monkeypatch.setattr(jobrecord, "cancel", lambda folder: signalled.append(folder) or True)
    _send_action(document, f"run-cancel:{run_id}:1")
    assert signalled == [str(root / "dep_fe607")]
    assert "Ended dep_fe607" in _footer(document)


def test_jobrecord_cancel_signals_only_a_verified_process(tmp_path, monkeypatch):
    """A recorded pid is signalled only when the process table shows it still running
    the recorded input; a Slurm id goes through scancel; the record then reads cancelled."""
    from tumult.synthesis import jobrecord, scheduler
    import socket
    folder = tmp_path / "run"
    folder.mkdir()
    (folder / "job.nml").write_text(
        "&job_record\n state='running', pid=%d, host='%s', input='x.nml',\n /\n"
        % (os.getpid(), socket.gethostname()))
    sent = []
    monkeypatch.setattr(jobrecord, "process_runs", lambda pid, marker: False)
    monkeypatch.setattr(jobrecord.os, "kill", lambda pid, sig: sent.append((pid, sig)))
    assert not jobrecord.cancel(folder) and not [s for s in sent if s[1] == 15]
    monkeypatch.setattr(jobrecord, "process_runs", lambda pid, marker: marker == "x.nml")
    assert jobrecord.cancel(folder)
    assert [s for s in sent if s[1] == 15] == [(os.getpid(), 15)]
    assert jobrecord.state_of(folder) == "cancelled"

    slurm = tmp_path / "slurm"
    slurm.mkdir()
    (slurm / "job.nml").write_text(
        "&job_record\n state='running', pid=1, host='node7', input='y.nml', slurm_job_id='99',\n /\n")
    cancelled = []
    monkeypatch.setattr(scheduler, "state", lambda job_id: "running")
    monkeypatch.setattr(scheduler, "scancel", lambda job_id: cancelled.append(job_id) or True)
    assert jobrecord.cancel(slurm) and cancelled == ["99"]
    assert jobrecord.state_of(slurm) == "cancelled"


def test_a_recovered_completion_is_imported_once_per_document(tmp_path, monkeypatch):
    """A reattached workflow delivers its result to the new document once; the record
    is then neither offered nor importable again there, while a third document that
    never saw it can import it explicitly, also once."""
    started, release, _ = _blocking_execute(monkeypatch)
    old_document, old_ctx, _ = _reconnecting_document(tmp_path, monkeypatch)
    old_delivered = []
    handle = _submit_synthesis(old_ctx, "t5750", _stub_request(tmp_path), old_delivered.append,
                                     meta={"type": "synth", "params": [["RT mode", "LTE"]]})
    assert started.wait(10)
    _destroy(old_document)

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    assert ctx.attach_workflow(handle.id) is not None
    release.set()
    assert ticks.settle(lambda: not ctx.jobs.active)
    assert [o.name for o in ctx.canvas.objects] == ["t5750"]
    assert ctx.canvas.objects[0].meta == {"type": "synth", "params": [["RT mode", "LTE"]]}
    assert ctx.workflow_refs.imported(handle.id) and not old_delivered
    assert ctx.discover_workflows() == [] and not len(ctx.recovery)
    assert ctx.import_workflow(handle.id) is None
    assert "already imported" in _footer(document)
    assert [o.name for o in ctx.canvas.objects] == ["t5750"]

    other_document, other, other_ticks = _reconnecting_document(tmp_path, monkeypatch)
    assert other.discover_workflows() == [handle.id]
    assert 'data-action="workflow-import"' in _recovery(other_document)
    assert "1 record to recover" in _processes(other_document)
    _send_action(other_document, f"workflow-import:{handle.id}:1")
    assert other_ticks.settle(lambda: not other.jobs.active)
    assert [o.name for o in other.canvas.objects] == ["t5750"]
    assert other.import_workflow(handle.id) is None
    assert [o.name for o in other.canvas.objects] == ["t5750"]


def test_slot_references_reattach_a_running_workflow_after_a_load(tmp_path, monkeypatch):
    """A save made while a workflow runs references it beside the workspace; loading that
    slot in a new document follows the workflow again under the new generation, and in
    the submitting document re-subscribes the existing row rather than duplicating it."""
    import store
    from tumult.spectrum_io import load_file

    started, release, _ = _blocking_execute(monkeypatch)
    old_document, old_ctx, old_ticks = _reconnecting_document(tmp_path, monkeypatch)
    path = _spectrum_file(tmp_path)
    old_ctx.canvas.add_object(load_file(path), name="obs", source=("file", str(path)))
    handle = _submit_synthesis(old_ctx, "saved", _stub_request(tmp_path), lambda r: None)
    assert started.wait(10)
    _session_action(old_document, "save with-workflow 1")
    slot = store.get_slot(1)
    assert slot["workflows"] == [{"id": handle.id, "folder": str(handle.folder),
                                  "name": "saved", "imported": False}]
    assert "workflows" not in slot["snapshot"]

    # the same document: the load boundary keeps the detached row and re-subscribes it
    report = old_ctx.restore_workspace(slot["snapshot"], "with-workflow", slot["workflows"])
    assert report.workflows == ["saved (running, followed)"]
    job, = old_ctx.jobs.active
    assert job.workflow_id == handle.id and job.generation is old_ctx.generation
    _destroy(old_document)

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    _session_load(document, 1, 2)
    assert "reattached" in _footer(document)
    job, = ctx.jobs.active
    assert job.workflow_id == handle.id and job.generation is ctx.generation
    release.set()
    assert ticks.settle(lambda: not ctx.jobs.active)
    assert [o.name for o in ctx.canvas.objects] == ["obs", "saved"]
    assert ctx.workflow_refs.imported(handle.id)
    _session_action(document, "save after 3")
    assert store.get_slot(2)["workflows"][0]["imported"] is True

    # the old document's row settled too, but delivered nothing to a closed document
    assert old_ticks.settle(lambda: not old_ctx.jobs.active)
    assert [o.name for o in old_ctx.canvas.objects] == ["obs"]


def _record_fixture(tmp_path, ident, state, **extra):
    """An on-disk workflow record as another server process would have left it."""
    import json
    folder = tmp_path / "synth" / "_workflows" / ident
    folder.mkdir(parents=True)
    log = tmp_path / "synth" / f"run-{ident}.log"
    log.write_text("dispatch log\n")
    run = tmp_path / "synth" / f"run-{ident}"
    record = {"version": 1, "id": ident, "state": state, "phase": state,
              "submitted": 1_700_000_000.0, "started": 1_700_000_001.0,
              "task": {"kind": "synthesis.synthesize", "args": {}}, "runs": [str(run)],
              "result": None, "label": f"{state}-job",
              "meta": None, **extra}
    (folder / "workflow.json").write_text(json.dumps(record))
    return folder


def test_on_disk_records_are_shown_by_their_true_state_after_a_restart(tmp_path, monkeypatch):
    """After a server restart, a record left running is interrupted (its logs are offered,
    nothing is resumed or signalled), a failed one keeps its error, and a finished one
    written by a previous server imports into this document."""
    from tumult.synthesis import submission, workflows

    interrupted = _record_fixture(tmp_path, "20260901T000000-interrupted", "running")
    _record_fixture(tmp_path, "20260901T000001-failed", "failed",
                    error="RuntimeError: dispatch exited with status 3", ended=1_700_000_100.0)
    _stub_synthesize(monkeypatch, lambda request, background_lines=True, cancel_token=None:
                     _stub_result(request))
    finished = submission.submit_operation("synthesis.synthesize", _stub_request(tmp_path),
                                           label="done-job")
    finished.result(timeout=10)

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    assert not len(ctx.recovery)                   # records wait for Recover
    _send_action(document, "workflow-discover:-:1")
    states = {entry.id: entry.state for entry in ctx.recovery}
    assert states == {"20260901T000000-interrupted": "interrupted",
                      "20260901T000001-failed": "failed", finished.id: "finished"}
    body = _recovery(document)
    assert "the process that ran it is gone" in body and "status 3" in body
    assert f"run-{interrupted.name}.log" in body
    assert body.count('data-action="workflow-import"') == 1
    assert 'data-action="workflow-attach"' not in body

    assert ctx.attach_workflow(interrupted.name) is None
    assert "cannot be attached" in _footer(document)
    assert submission.open_workflow(interrupted).cancel() is False
    assert ctx.import_workflow(interrupted.name) is not None
    assert ticks.settle(lambda: ctx.jobs.retained)
    assert "is interrupted" in ctx.jobs.retained[0].error and not ctx.canvas.objects

    _send_action(document, f"workflow-import:{finished.id}:2")
    assert ticks.settle(lambda: ctx.canvas.objects)
    assert [o.name for o in ctx.canvas.objects] == ["done-job"]
    _send_action(document, "workflow-forget:20260901T000001-failed:3")
    assert [entry.id for entry in ctx.recovery] == []
    assert "status 3" not in _recovery(document)


def test_a_slot_whose_workflow_record_is_gone_reports_it(tmp_path, monkeypatch):
    """A referenced workflow whose folder no longer exists is reported as skipped, with
    its folder, and offered nowhere; the reference itself is kept for the next save."""
    import store

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    _with_spectrum(ctx)
    _session_action(document, "save plain 1")
    gone = store.workflow_reference("20260901T000000-gone", str(tmp_path / "nowhere"), "gone")
    store.save_session("dangling", store.get_session(1), [gone])
    _session_load(document, 2, 2)
    footer = _footer(document)
    assert "workflow gone" in footer and "no record" in footer
    assert str(tmp_path / "nowhere") in footer
    assert not len(ctx.recovery) and "20260901T000000-gone" in ctx.workflow_refs


def _fire_timeouts(document):
    """Run the document's pending timeout callbacks (the debounced autosave) now."""
    for callback in list(document.session_callbacks):
        # an earlier callback of this pass may have replaced a later one (a debounce)
        if type(callback).__name__ == "TimeoutCallback" and callback in document.session_callbacks:
            document.remove_timeout_callback(callback)
            callback.callback()


def test_the_workspace_is_autosaved_after_a_change_and_on_close(tmp_path, monkeypatch):
    """A fresh document writes no autosave (a blank workspace must not replace the last
    real one); the first change writes one slot, later changes update that slot in place,
    it lists first as the autosave, and closing the tab writes the last state at once."""
    import store

    document, ctx = _build_document()
    _fire_timeouts(document)
    assert store.autosave_slot() is None

    from tumult.spectrum_io import load_file
    obs_path = _spectrum_file(tmp_path)
    ctx.canvas.add_object(load_file(obs_path), name="obs", source=("file", str(obs_path)))
    _fire_timeouts(document)
    slot = store.autosave_slot()
    assert slot is not None and slot["name"] == store.AUTOSAVE_NAME
    assert len(slot["snapshot"]["spectra"]) == 1
    sid = slot["id"]

    _session_action(document, "save named 1")                # a manual save is another slot
    x = np.arange(5000., 5010., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="second", source=("synth", None))
    _fire_timeouts(document)
    slots = store.list_sessions()
    assert [s["id"] for s in slots][0] == sid and slots[0]["autosave"]
    assert len(store.autosave_slot()["snapshot"]["spectra"]) == 2
    assert len(slots) == 2

    ctx.canvas.remove_object(ctx.canvas.objects[-1].uid)
    for callback in list(document.session_destroyed_callbacks):
        callback(None)                                       # the tab closes: no debounce
    assert len(store.autosave_slot()["snapshot"]["spectra"]) == 1

    # the next document resumes from it
    document, ctx = _build_document()
    _session_load(document, sid, 3)
    assert [o.name for o in ctx.canvas.objects] == ["obs"]


# --- synthesis and fit requests collected by the modal -----------------------------------
_GRID = {"id": "g1", "name": "Test grid", "folder": "", "pattern": "", "format": "Marcs"}


def _synth_document(tmp_path, monkeypatch, run=None, grid_id="g1", executed=None):
    """A document whose synthesis setup resolves without files on disk: one registered
    grid (or the Individual pseudo-grid), no line list, and the workflow root under
    ``tmp_path``. ``run`` overrides the ``synth_run`` section; ``executed`` collects the
    resolved request of every workflow the stubbed execution receives."""
    import store
    store.set_section("synthesis", {
        "dispatch_cwd": str(tmp_path), "synth_output_dir": "synth", "linelist": "",
        "abund_file": "", "absmet": "", "atmos_grids": [dict(_GRID)],
        "atmos_grid_id": grid_id, "ref_grid_id": "g1"})
    store.set_section("synth_run", {"teff": 5750.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0,
                                    "mode": "LTE", "abundances": "", "atoms": [],
                                    "relative_abundances": True, **(run or {})})
    if executed is not None:
        def execute(request, background_lines=True, cancel_token=None):
            executed.append(request)
            return _stub_result(request)
        _stub_synthesize(monkeypatch, execute)
    document, ctx = _build_document()
    messages = []
    ctx.notify = lambda text, level="info", **_: messages.append(text)
    return document, ctx, messages


def _departure_folder(tmp_path):
    folder = tmp_path / "dep"
    (folder / "save").mkdir(parents=True)
    (folder / "save" / "pop_Mg_meta.txt").write_text("")
    (folder / "params.nml").write_text(
        "&atmos_params\n atmos_file='one.mod' nx=1 ny=1 nz=6\n/\n"
        "&atom_params\n atom_file='/atoms/mg.atom'\n/\n"
        "&composition_params\n abundances=''\n/\n")
    return folder


_presses = itertools.count()


def _submit_modal(document, ctx, ticks):
    """Press the modal's footer button through the document's action channel and run the
    next tick — the actual handoff to the daemon is deferred there now, so the dialog's
    close patch can reach the browser before that (possibly slow) round trip, not after
    it. Returns the number of jobs the document admitted for it."""
    before = len(ctx.jobs)
    _send_action(document, f"synth-submit:{next(_presses)}")
    ticks.run()
    return len(ctx.jobs) - before


def test_free_range_synthesis_is_one_durable_window_workflow(tmp_path, monkeypatch):
    """A free range is one window, submitted as a detached
    durable workflow whose record carries the modal's name and Info-tab meta; its result
    reaches the canvas as a synth object."""
    executed = []
    document, ctx, messages = _synth_document(tmp_path, monkeypatch, executed=executed)
    ticks = _Ticks(document, monkeypatch)
    modal = ctx.synth_modal
    modal.open_range(5000.0, 5001.0, "synth")
    assert _submit_modal(document, ctx, ticks) == 1
    job, = ctx.jobs.active
    assert job.detached and job.handle is not None and job.name == "t5750g4.4m+0.0v1.0"
    assert ticks.settle(lambda: not ctx.jobs.active)
    request, = executed
    assert request.windows == ((5000., 5001.),)
    assert isinstance(request.atmosphere, AtmosphereGrid)
    assert request.star.tuple() == (5750., 4.4, 0., 1.)
    assert request.departures is None
    obj, = ctx.canvas.objects
    assert obj.name == "t5750g4.4m+0.0v1.0" and obj.source == ("synth", None)
    assert ["Window", "5000.00–5001.00 Å"] in obj.meta["params"]
    assert job.handle.meta == obj.meta and job.handle.label == obj.name
    assert not messages


def test_mask_set_synthesis_keeps_fit_and_synth_windows_apart(tmp_path, monkeypatch):
    executed = []
    document, ctx, messages = _synth_document(tmp_path, monkeypatch, executed=executed)
    ticks = _Ticks(document, monkeypatch)
    x = np.arange(5000., 5010., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="obs",
                          source=("file", "/nowhere/obs.fits"))
    range_tool = ctx.tools._tools["range"]
    ctx.tools.toggle("range")
    range_tool._on_select(5004.5, 5005.5)
    range_tool.new_mask()
    mset = ctx.masks.mset
    mset.masks[0].include = True
    ctx.synth_modal.open_masks(mset, "synth")
    assert _submit_modal(document, ctx, ticks) == 1
    assert ticks.settle(lambda: not ctx.jobs.active)
    request, = executed
    assert request.windows == ((5004.0, 5006.0),)
    synth, = [o for o in ctx.canvas.objects if o.source[0] == "synth"]
    assert ["Line mask", mset.name] in synth.meta["params"]


def test_fixed_file_snapshot_and_nlte_requests_are_captured_from_the_setup(tmp_path,
                                                                            monkeypatch):
    """The same context() yields a picked 1D file, a 3D snapshot of the Individual
    pseudo-grid and an NLTE request whose saved departures is forwarded; each is
    submitted as a durable workflow and the field survives the record's JSON round trip."""
    from tumult.synthesis import Atmosphere, Snapshot, run_values
    one = tmp_path / "one.mod"
    one.write_text("")
    executed = []
    document, ctx, messages = _synth_document(
        tmp_path, monkeypatch, executed=executed,
        run={"atmos_mode": "pick", "atmos_file": str(one)})
    ticks = _Ticks(document, monkeypatch)
    modal = ctx.synth_modal
    modal.open_range(5000.0, 5001.0, "synth")
    assert _submit_modal(document, ctx, ticks) == 1
    assert ticks.settle(lambda: not ctx.jobs.active)
    picked = executed[-1]
    assert isinstance(picked.atmosphere, Atmosphere)
    assert picked.atmosphere.path == one and picked.atmosphere.format == "Marcs"

    import store
    store.update_pref("synthesis", "atmos_grid_id", "__individual__")
    store.update_pref("synth_run", "atmos_individual_format", "stagger")
    modal.setup.discard()        # re-read the store: the setup holds a draft of it
    assert modal.setup.is_snapshot
    modal.open_range(5000.0, 5001.0, "synth")
    assert _submit_modal(document, ctx, ticks) == 1
    assert ticks.settle(lambda: not ctx.jobs.active)
    cube = executed[-1]
    assert isinstance(cube.atmosphere, Snapshot) and cube.atmosphere.format == "stagger"
    assert cube.atmosphere.resolution == (30, 30, 128) and cube.windows == ((5000., 5001.),)

    store.update_pref("synth_run", "atmos_individual_format", "Marcs")
    modal.setup.discard()
    modal.mode_form.set_values(mode="NLTE", departures=str(_departure_folder(tmp_path)))
    modal.atoms.selected = ["/atoms/mg.atom"]
    modal.open_range(5000.0, 5001.0, "synth")
    assert _submit_modal(document, ctx, ticks) == 1
    job, = ctx.jobs.active
    assert ticks.settle(lambda: not ctx.jobs.active)
    nlte = executed[-1]
    assert run_values(nlte)["mode"] == "NLTE" and run_values(nlte)["atoms"] == ["/atoms/mg.atom"]
    assert nlte.departures is not None and nlte.departures.species == ("Mg",)
    assert ["Saved departures", str(nlte.departures.folder)] in job.handle.meta["params"]
    assert ["Model atoms", "mg.atom"] in job.handle.meta["params"]
    assert [o.name for o in ctx.canvas.objects][-1] == "t5750g4.4m+0.0v1.0·NLTE"
    assert not messages


def test_unsupported_requests_are_refused_before_submission(tmp_path, monkeypatch):
    """A fit of a 3D snapshot is refused with the reason and admits no job. A saved
    saved departures are not offered outside NLTE synthesis (the field is hidden for LTE and
    for fits), so a stored one is ignored there rather than refused."""
    from tumult.synthesis import operations
    one = tmp_path / "one.mod"
    one.write_text("")
    document, ctx, messages = _synth_document(
        tmp_path, monkeypatch, grid_id="__individual__",
        run={"atmos_file": str(one), "departures": str(_departure_folder(tmp_path)),
             "atmos_individual_format": "stagger"})
    monkeypatch.setattr(operations, "synthesize",
                        lambda *a, **k: pytest.fail("a refused request was executed"))
    monkeypatch.setattr(operations, "fit", lambda *a, **k: pytest.fail("a refused fit ran"))
    ticks = _Ticks(document, monkeypatch)
    modal = ctx.synth_modal
    x = np.arange(4999., 5002., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="obs",
                          source=("file", "/nowhere/obs.fits"))
    modal.fit.target = ctx.canvas.objects[0].uid
    modal.fit.abund_rows = {"Mg": 0.0}

    modal.open_range(5000.0, 5001.0, "fit")          # a fit of a 3D snapshot
    assert modal.setup.mode_form.fields["departures"].hidden
    assert _submit_modal(document, ctx, ticks) == 0
    assert messages[-1] == ("Cannot fit: Fitting is available for 1D atmospheres and "
                            "grids only.")
    assert not len(ctx.jobs)
    # the dialog stayed open (nothing re-renders its footer to un-stick the submit
    # button on its own), so the reject pulse is what resets it client-side
    assert modal.dialog.value.startswith("reject:")


def test_fit_is_a_durable_workflow_built_from_the_same_request(tmp_path, monkeypatch):
    """A fit gives the setup's request the job's windows and runs as a detached
    ``synthesis.fit`` record (it survives the session closing, its record names the
    request and the cropped observation) and adds the raw best-fit spectrum with the
    broadening as its live convolution and the request's provenance."""
    import threading
    from tumult.synthesis import operations
    from tumult.synthesis.request import SpectrumFit
    document, ctx, messages = _synth_document(tmp_path, monkeypatch)
    ticks = _Ticks(document, monkeypatch)
    modal = ctx.synth_modal
    x = np.arange(4999., 5003., .01)
    ctx.canvas.add_object(Spectrum(x, np.ones_like(x)), name="obs",
                          source=("file", "/nowhere/obs.fits"))
    modal.fit.target = ctx.canvas.objects[0].uid
    modal.fit.abund_rows = {"Mg": 0.2}
    modal.fit.flags["vshift"] = True
    modal.fit.flags["teff"] = True
    fitted, started, release = [], threading.Event(), threading.Event()

    def fit(request, *, target, cancel_token=None, **kwargs):
        fitted.append((target, request, kwargs, cancel_token))
        started.set()
        assert release.wait(10)
        cancel_token.check()
        xx = np.linspace(5000., 5001., 5)
        return SpectrumFit({"teff": 5600., "Mg": 0.35, "vshift": 1.5},
                           Spectrum(xx, np.ones_like(xx)), {"vshift": 1.5}, request)

    monkeypatch.setattr(operations, "fit", fit)
    modal.open_range(5000.0, 5001.0, "fit")
    assert _submit_modal(document, ctx, ticks) == 1
    job, = ctx.jobs.active
    assert job.name == "fit·obs" and job.heavy and job.detached
    assert job.handle is not None and job.handle.kind.name == "synthesis.fit"
    assert started.wait(10)
    observation, request, kwargs, token = fitted[0]
    assert kwargs["fit_stellar"] == ["teff"] and kwargs["abund_fits"] == {"Mg": 0.2}
    assert kwargs["fit_post"] == {"vshift": 0.0}
    assert request.windows == ((5000., 5001.),) and kwargs["fit_windows"] == [[5000., 5001.]]
    assert isinstance(request.atmosphere, AtmosphereGrid)
    assert request.config.root == tmp_path.resolve()
    ox, _ = observation.values                  # the record carries the fitted span only
    assert ox.min() >= 5000. and ox.max() <= 5001. and len(ox) > 50
    release.set()
    assert ticks.settle(lambda: not ctx.jobs.active)
    synth, = [o for o in ctx.canvas.objects if o.source[0] == "synth"]
    assert synth.name == "fit·obs" and synth.conv == {"vshift": 1.5}
    assert synth.meta["stellar"]["teff"] == 5600.
    assert ["Atmosphere", "Test grid grid, interpolated"] in synth.meta["params"]
    assert ["Abundances ([X/Fe])", "Mg=+0.35"] in synth.meta["params"]
    assert modal.run_form.value("teff") == 5600.        # the best fit becomes the setup
    assert job.handle.id in ctx.workflow_refs and ctx.workflow_refs.imported(job.handle.id)

    # a second fit keeps running when the session closes; its record delivers the result
    started.clear()
    release.clear()
    modal.open_range(5000.0, 5001.0, "fit")
    assert _submit_modal(document, ctx, ticks) == 1
    job, = ctx.jobs.active
    assert started.wait(10)
    for callback in list(document.session_destroyed_callbacks):
        callback(None)
    assert not fitted[1][3].cancelled
    release.set()
    result = job.handle.result(timeout=10)
    assert result.best["teff"] == 5600. and result.post == {"vshift": 1.5}


# --- the Processes panel's rows and the job details modal ----------------------------------
class _Live:
    """Stands in for ``daemon.LiveJobs``: the daemon's live jobs and its thread use."""

    def __init__(self, jobs, used=0, budget=8):
        self._jobs, self.used, self.budget = {j["job"]: j for j in jobs}, used, budget

    def jobs(self):
        return list(self._jobs.values())

    def get(self, job):
        return self._jobs.get(job)

    def external_threads(self):
        return sum(j.get("threads") or 0 for j in self._jobs.values() if j["type"] == "external")


def test_panel_rows_name_the_kind_and_how_each_job_runs(monkeypatch):
    """Two lines per row: ``Kind · name``, then threads, queued, stopping… or failed; the
    daemon's jobs of other sessions and external runs follow with details and cancel,
    the header counts external threads over the budget, and a job admitted a moment ago
    is not drawn."""
    import time
    import processes
    from jobs import Job, JobTable
    monkeypatch.setattr(processes, "SHOW_AFTER_S", 0.5)
    now = time.time()
    table = JobTable()
    strengths = table.add(Job(1, "line-strengths", CancelToken(), None, submitted=now - 60,
                              started=now - 30, state="running"))
    table.add(Job(2, "obs.txt", CancelToken(), None, path="/data/obs.txt", submitted=now))
    table.add(Job(3, "abund-synth:Fe", CancelToken(), None, submitted=now - 5,
                  state="failed", error="RuntimeError: no lines"))
    live = _Live([
        {"job": "w1", "type": "workflow", "kind": "synthesis.synthesize", "label": "Sun 3D NLTE",
         "state": "running", "threads": 8, "started": now - 720, "folder": "/r/w1"},
        {"job": "w2", "type": "workflow", "kind": "sparams.ionbalance", "label": "ionbal:Ti",
         "state": "queued", "folder": "/r/w2"},
        {"job": "run-abc", "type": "external", "kind": "dispatch", "label": "dep_fe",
         "state": "running", "threads": 3, "started": now - 20, "folder": "/r/dep_fe"}],
        used=8, budget=16)
    body = processes.processes_html(table, live, now=now)
    assert "Line strengths" in body and "1 thread<" in body
    assert f'data-since="{strengths.started:.3f}"' in body
    assert "obs.txt" not in body                        # admitted just now: not drawn
    assert "Abundance synthesis · Fe" in body and "failed · RuntimeError: no lines" in body
    assert "Synthesis · Sun 3D NLTE" in body and "8 threads" in body
    assert "Ionisation balance · Ti" in body and ">queued<" in body
    assert "Dispatch run · dep_fe" in body and 'data-action="run-cancel" data-uid="run-abc"' in body
    assert body.count('data-action="job-details"') == 4    # failed, w1, w2, the external run
    assert processes.usage_html(live) == "11/16 threads"
    assert "Import · obs.txt" == processes.job_title("import obs.txt")
    assert "Load · obs.txt" in processes.processes_html(table, None, now=now + 1)


def test_the_details_modal_shows_a_workflow_record(tmp_path, monkeypatch):
    """Progress, execution, timing, the submitter's inputs, a log tail and the actions a
    live workflow of another session allows."""
    import json
    import processes
    from tumult.synthesis import submission
    folder = tmp_path / "synth" / "_workflows" / "w1"
    folder.mkdir(parents=True)
    run = tmp_path / "synth" / "run_a"
    run.with_suffix(".log").write_text("".join(f"line {i}\n" for i in range(40)))
    (folder / "workflow.json").write_text(json.dumps({
        "id": "w1", "state": "queued", "phase": "queued", "submitted": 1_700_000_000.0,
        "label": "t5750", "runs": [str(run)], "task": {"kind": "synthesis.synthesize", "args": {}},
        "meta": {"type": "synth", "params": [["RT mode", "LTE"]]},
        "driver": {"kind": "daemon", "host": None, "pid": None}}))
    handle = submission.open_workflow(folder)
    body = processes.workflow_details_html(handle, {"threads": 0, "kind": "synthesis.synthesize"})
    assert "Synthesis · t5750" in body and "RT mode" in body and "LTE" in body
    assert "daemon" in body and "2023-11-14" in body
    assert "line 39" in body and "line 19" not in body       # the tail only
    assert 'data-action="workflow-attach"' in body and 'data-action="workflow-cancel"' in body
    assert 'data-copy="w1"' in body


@pytest.mark.parametrize('tool_name', ['abundance', 'ionbal'])
def test_measurement_continuum_markers_keep_target_coordinates(tool_name):
    """The real document wires each measurement's checkbox and vertical edge markers."""
    from tumult.abundance import AbundanceResult, LineResult
    from tumult.measurement import BalanceResult
    from tumult.linemask import LineMask
    document, ctx = _build_document()
    xx = np.linspace(4999., 5001., 401)
    yy = (1.03 + .015 * (xx - 5000)) * (1 - .3*np.exp(-.5*((xx-5000)/.1)**2))
    obj = ctx.canvas.add_object(Spectrum(xx, yy), name='target')
    # The prepared target has a deliberately non-flat global continuum.
    obj.norm_knots = (np.array([4999., 5001.]), np.array([1., 1.01]))
    ctx.canvas.normalise_object(obj)
    before = obj.spectrum.ydata.copy()
    knots = tuple(v.copy() for v in obj.norm_knots)
    ln = LineResult(5000., 'Fe I', 1., -1.,
                    mask=LineMask(lmin=4999.5, lmax=5000.5, smin=4999., smax=5001.))
    cls = AbundanceResult if tool_name == 'abundance' else BalanceResult
    extra = {'symbol': 'Fe', 'summary': {}} if tool_name == 'abundance' else {}
    res = cls(lines=[ln], element='Fe', ref_uid=obj.uid, profile='gauss',
               snr=100., cont_pct=.5, **extra)
    tool = ctx.tools._tools[tool_name]
    tool.results.append(res)
    ctx.tools.activate(tool_name)
    tool.select_result(res.uid)
    tool.actions.value = 'continuum:0:1'
    assert ln.renorm and tool.continuum_markers.visible
    assert 'line-continuum-chk' in tool.table_html(res)
    ew = ln.ew
    # An attempted horizontal drag is snapped back to the synth-window edge.
    tool.continuum_source.data = dict(x=[4999.2, 5000.8],
                                      y=[p[1] + .01 for p in ln.continuum])
    assert [p[0] for p in ln.continuum] == [4999., 5001.]
    assert ln.ew > ew
    np.testing.assert_array_equal(obj.spectrum.ydata, before)
    for original, current in zip(knots, obj.norm_knots):
        np.testing.assert_array_equal(original, current)
    tool.deactivate()
    assert not tool.continuum_markers.visible
    for callback in list(document.session_destroyed_callbacks):
        callback(None)
