"""A running analysis result: listed and in view at once, a progress strip that says how far
it is, and its first row selected when the lines arrive.

What matters: someone who submits an abundance run (the balance and correction views
share the same base) sees the new object selected and its view opened before the
calculation has done anything; while it runs, the table pane shows the stage under way
and the steps done over those expected (with the elapsed time and an estimate of what is
left); and when a snapshot with lines lands, the first row of the table as displayed is
selected, so the line panels show something rather than nothing. The expected answers
follow from the tool's own contract: the submission is captured instead of reaching a
daemon, and the snapshots are hand-built results as :func:`tumult.abundance.fit` reports
them. A defect this catches: the object only appearing when the measurement returns, a
strip that outlives the job, or a selected result opening on no row.
"""
import re

import numpy as np
import pytest

from tests.support.spectra import gaussian_line

pytest.importorskip("bokeh")


class _Submission:
    """Stands in for ``ctx.submit_operation``: keeps the callbacks a run hands over."""

    def __init__(self):
        self.calls = []

    def __call__(self, name, operation, request, on_done, *, kwargs=None, meta=None,
                 on_progress=None, on_error=None, on_native=None):
        from tumult.synthesis.cancel import CancelToken

        class Job:
            token = CancelToken()
        self.calls.append(dict(name=name, operation=operation, request=request,
                               kwargs=kwargs, done=on_done, progress=on_progress,
                               native=on_native))
        return Job()


def _line(lam0, a=float("nan")):
    from tumult.abundance import LineResult
    from tumult.linemask import LineMask
    return LineResult(lam0=lam0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=50.0, a=a,
                      mask=LineMask(lam0=lam0, lmin=lam0 - 0.2, lmax=lam0 + 0.2,
                                    smin=lam0 - 0.5, smax=lam0 + 0.5))


def _snapshot(res, lines, stage, done, total):
    from copy import deepcopy
    from tumult.measured import report_progress
    snap = deepcopy(res)
    snap.lines = lines
    report_progress(snap, stage, done, total)
    return snap


@pytest.fixture
def studio(bokeh_document, studio_config):
    from TumultStudio import app
    return app.build()


def test_run_lists_the_object_at_once_and_streams_its_progress(studio, monkeypatch):
    from tumult.abundance import AbundanceResult, aggregate
    from tumult.spectrum import Spectrum

    ctx = studio
    tool = ctx.tools.get("abundance")
    submission = _Submission()
    monkeypatch.setattr(ctx, "submit_operation", submission)
    monkeypatch.setattr(ctx, "schedule", lambda fn: fn())     # the deferred submission
    x, y = gaussian_line()
    spec = Spectrum(x + 5000.0, y)
    monkeypatch.setattr(tool, "get_observation", lambda uid: spec)

    res = AbundanceResult(element="Fe", symbol="Fe", lines=[], summary=aggregate([]),
                          stage="lines", uid="abc12345", ref_uid="obs-1", method="cog")
    tool._add_result(res, "Fe abundance")
    tool.select_result(res.uid)
    tool.run_synth(object(), spec, (100.0, 0.5), line_mask=object())

    # the object is listed, selected and in view before anything has been computed
    assert tool.selected() is res and tool.is_focused()
    assert ctx.workspace.focus_kind == tool.kind
    row = tool.results_html()
    assert 'data-result="abc12345"' in row and "Measuring windows · 0%" in row
    strip = tool._progress_bridge.value
    assert "result-progress" in strip and "Measuring windows" in strip
    assert "No Fe lines measured yet" in tool._table_bridge.value
    assert tool.jumped is None
    # the submission continues that object and carries the mask to measure into it
    [call] = submission.calls
    assert call["operation"] == "abundance.fit"
    assert call["kwargs"]["measured"] is res and call["kwargs"]["line_mask"] is not None

    # the measured lines arrive: the first row of the table (sorted by wavelength, as
    # displayed) is selected, and the strip says where the run is
    call["progress"](_snapshot(res, [_line(5001.0), _line(5000.0)], "Screening lines", 1, 5))
    table = tool._table_bridge.value
    assert len(res.lines) == 2 and tool.jumped == 1          # the 5000 Å line
    rows = re.findall(r'<tr class="abund-row([^"]*)" data-line="(\d)"', table)
    assert rows == [(" row-selected", "1"), ("", "0")]
    strip = tool._progress_bridge.value
    assert re.search(r'data-done="1" data-total="5"', strip)
    assert "Screening lines" in strip and "1/5" in strip and "left" in strip
    assert "Screening lines · 20%" in tool.results_html()

    # the running synthesis' own fraction (its RT patch tasks) fills the step under way
    # without repainting the table
    table = tool._table_bridge.value
    call["native"](0.5)
    strip = tool._progress_bridge.value
    assert 'data-run="0.50"' in strip and "synthesis 50%" in strip
    assert 'style="width: 30.0%"' in strip and tool._table_bridge.value is table
    assert "Screening lines · 30%" in tool.results_html()
    # what is shown never moves back: a run's estimate slipping, the next run of the
    # same step starting from nothing, or a snapshot whose expected total grew all hold
    # the bar where it was until the work catches up
    call["native"](0.2)
    assert 'data-run="0.20"' in tool._progress_bridge.value
    assert 'style="width: 30.0%"' in tool._progress_bridge.value
    call["native"](None)
    assert "data-run" not in tool._progress_bridge.value
    assert 'style="width: 30.0%"' in tool._progress_bridge.value
    call["progress"](_snapshot(res, [_line(5001.0), _line(5000.0)], "Screening lines", 1, 8))
    strip = tool._progress_bridge.value
    assert 'data-total="8"' in strip and 'style="width: 30.0%"' in strip
    call["progress"](_snapshot(res, [_line(5001.0, 7.4), _line(5000.0, 7.5)],
                               "Curve of growth · round 2", 3, 5))
    strip = tool._progress_bridge.value
    assert 'data-done="3"' in strip and 'style="width: 60.0%"' in strip
    assert "data-run" not in strip

    # a later round kept the selection; the finished result drops the strip
    assert tool.jumped == 1
    final = _snapshot(res, [_line(5001.0, 7.4), _line(5000.0, 7.5)], "", 5, 5)
    final.progress, final.stage = {}, "uncert"
    call["done"](final)
    assert tool.selected() is res and res.stage == "uncert"
    assert tool._progress_bridge.value == ""
    assert "%" not in tool.results_html()


def test_a_selected_result_opens_on_its_first_row(studio):
    from tumult.abundance import AbundanceResult, aggregate

    tool = studio.tools.get("abundance")
    res = AbundanceResult(element="Fe", symbol="Fe", lines=[_line(5002.0), _line(5001.0)],
                          summary=aggregate([]), stage="uncert", uid="def67890",
                          method="cog")
    tool._add_result(res, "Fe abundance")
    tool.select_result(res.uid)
    assert tool.jumped == 1
    assert 'class="abund-row row-selected" data-line="1"' in tool._table_bridge.value
    # through the jump itself, which moves the view onto that line's window
    assert tool._last_span_aa == pytest.approx(1.0)      # its synth window
    # selecting again opens on the first row again; a click elsewhere is kept until then
    tool.jump_line("0")
    assert tool.jumped == 0
    tool.select_result(res.uid)
    assert tool.jumped == 1


def test_a_followed_job_polls_its_running_synthesis(studio, monkeypatch):
    """While a job runs, the document reads how far its dispatch runs are every poll and
    hands each change to the submitter; the poll ends with the job. The handle is a
    stand-in whose fraction the test moves; the document's periodic callbacks and the
    IO pool are driven here instead of by a server."""
    import time
    from TumultStudio.execution import workflows as wf
    from TumultStudio.execution.jobs import FINISHED, Job, RUNNING
    from tumult.synthesis.cancel import CancelToken

    ctx = studio
    service = ctx.workflows
    monkeypatch.setattr(ctx, "schedule", lambda fn: fn())

    class Handle:
        fraction = 0.25
        kind = type("Kind", (), {"progress": False})()
        detached = "daemon"

        def native_fraction(self):
            return self.fraction

    class Pool:
        def submit(self, fn):
            fn()
    monkeypatch.setattr(wf, "IO_POOL", Pool())

    seen = []
    job = Job(pid=1, name="run", token=CancelToken(), generation=ctx.generation, state=RUNNING)
    service._follow_progress(job, Handle(), None, seen.append)
    ticks = [cb for cb in ctx.doc.session_callbacks if getattr(cb, "period", None) == wf.NATIVE_POLL_MS]
    assert len(ticks) == 1
    tick = ticks[0].callback
    tick()
    tick()                                   # unchanged: not reported again
    assert seen == [0.25]
    Handle.fraction = 0.5
    tick()
    Handle.fraction = None                   # between runs
    tick()
    assert seen == [0.25, 0.5, None]
    job.state = FINISHED
    tick()
    assert seen == [0.25, 0.5, None]
    assert not [cb for cb in ctx.doc.session_callbacks if getattr(cb, "period", None) == wf.NATIVE_POLL_MS]
