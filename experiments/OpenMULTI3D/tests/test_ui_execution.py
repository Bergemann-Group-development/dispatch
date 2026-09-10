"""Background computations must not apply stale UI state or wait on HOP locks."""
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from bokeh.models import ColumnDataSource, Range1d, TapTool

from tumult import Spectrum
from tumult.synthesis.cancel import CancelToken


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
    tool.ctx = SimpleNamespace(submit=lambda name, run, done: queued.append((run, done)))

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


def test_result_jobs_snapshot_progress_and_reject_stale_completion():
    from jobs import ResultJobs
    queued, callbacks, rendered, completed = [], [], [], []
    result = SimpleNamespace(uid='fit', lines=[SimpleNamespace(a=0.)], termination='measured')
    results = {'fit': result}
    ctx = SimpleNamespace(submit=lambda name, run, done: queued.append((run, done)),
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
    ctx = SimpleNamespace(submit=lambda name, run, done: queued.append((run, done)),
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
