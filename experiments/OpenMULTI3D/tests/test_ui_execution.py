"""Background computations must not apply stale UI state or wait on HOP locks."""
import html
from pathlib import Path
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

    def refuse(request, **kwargs):
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
    _session_action(document, "load 1 2")
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

    _session_action(document, "load 1 3")
    state = _state(ctx)
    assert [o[1] for o in state["objects"]] == ["obs", "model"]
    assert [m[1] for m in state["masks"]] == ["fe-lines"]
    assert state["results"]["abundance"][0] == [("res-a", ids.obs)]
    assert state["focus"] == ("corr", "abundcorr")

    # the empty slot: every collection replaced by nothing, the plain canvas back
    store.save_session("empty", {"schema_version": 2})
    _session_action(document, "load 3 4")
    assert "Loaded “empty”: 0 spectra, 0 line masks, 0 atoms, 0 results" in _footer(document)
    state = _state(ctx)
    assert state["objects"] == [] and state["masks"] == [] and state["atoms"] == []
    assert all(r[0] == [] and r[1] is None for r in state["results"].values())
    assert state["focus"] == ("spectrum", None) and state["selected"] is None
    # and the other way round: loading the "other" slot into the emptied document
    _session_action(document, "load 2 5")
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
    _session_action(document, "load 2 2")
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
    _session_action(document, "load 3 3")
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
        _session_action(document, f"load {sid} {sid}")
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
    _session_action(document, "load 2 2")
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
    from tumult.synthesis import workflows
    started, release, tokens = threading.Event(), threading.Event(), []

    def execute(resolved, cancel_token=None):
        tokens.append(cancel_token)
        started.set()
        assert release.wait(10)
        cancel_token.check()
        return _stub_result(resolved)

    monkeypatch.setattr(workflows, "execute_request", execute)
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
    handle = old_ctx.submit_workflow("first", _stub_request(tmp_path), old_delivered.append)
    assert started.wait(10)
    assert old_ctx.workflow_refs.to_list() == [
        {"id": handle.id, "folder": str(handle.folder), "name": "first", "imported": False}]
    _destroy(old_document)

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    entry, = ctx.recovery              # offered at build: the server still runs it
    assert entry.id == handle.id and entry.state == "running" and entry.name == "first"
    assert 'data-action="workflow-attach"' in _processes(document)
    assert not ctx.jobs.active
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


def test_a_recovered_completion_is_imported_once_per_document(tmp_path, monkeypatch):
    """A reattached workflow delivers its result to the new document once; the record
    is then neither offered nor importable again there, while a third document that
    never saw it can import it explicitly, also once."""
    started, release, _ = _blocking_execute(monkeypatch)
    old_document, old_ctx, _ = _reconnecting_document(tmp_path, monkeypatch)
    old_delivered = []
    handle = old_ctx.submit_workflow("t5750", _stub_request(tmp_path), old_delivered.append,
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
    assert 'data-action="workflow-import"' in _processes(other_document)
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
    handle = old_ctx.submit_workflow("saved", _stub_request(tmp_path), lambda r: None)
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
    _session_action(document, "load 1 2")
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
              "request": {}, "runs": [str(run)], "result": None, "label": f"{state}-job",
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
    monkeypatch.setattr(workflows, "execute_request",
                        lambda resolved, cancel_token=None: _stub_result(resolved))
    finished = submission.submit(_stub_request(tmp_path), label="done-job")
    finished.result(timeout=10)
    assert finished.id not in submission._live

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    assert not len(ctx.recovery)                   # nothing is live in this server
    _send_action(document, "workflow-discover:-:1")
    states = {entry.id: entry.state for entry in ctx.recovery}
    assert states == {"20260901T000000-interrupted": "interrupted",
                      "20260901T000001-failed": "failed", finished.id: "finished"}
    body = _processes(document)
    assert "the server that ran it is gone" in body and "status 3" in body
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
    assert "status 3" not in _processes(document)


def test_a_slot_whose_workflow_record_is_gone_reports_it(tmp_path, monkeypatch):
    """A referenced workflow whose folder no longer exists is reported as skipped, with
    its folder, and offered nowhere; the reference itself is kept for the next save."""
    import store

    document, ctx, ticks = _reconnecting_document(tmp_path, monkeypatch)
    _with_spectrum(ctx)
    _session_action(document, "save plain 1")
    gone = store.workflow_reference("20260901T000000-gone", str(tmp_path / "nowhere"), "gone")
    store.save_session("dangling", store.get_session(1), [gone])
    _session_action(document, "load 2 2")
    footer = _footer(document)
    assert "workflow gone" in footer and "no record" in footer
    assert str(tmp_path / "nowhere") in footer
    assert not len(ctx.recovery) and "20260901T000000-gone" in ctx.workflow_refs
