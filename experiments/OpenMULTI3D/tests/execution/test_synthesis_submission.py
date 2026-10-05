"""Durable submission of a synthesis request through the workflow daemon.

The daemon runs in a thread of the test process (the ``workflow_daemon`` fixture), so
the stubbed synthesis is the one it executes.
"""
import numpy as np

from tumult import Spectrum
from tumult.synthesis import Atmosphere, SynthConfig, SynthesisResult, SynthRequest


def _request(atmosphere, tmp_path, windows=(), star=(5750, 4.5, 0, 1), **config):
    config = SynthConfig(**{'root': tmp_path, 'cache': 'out', 'R_sampling': 50000,
                            **config})
    return SynthRequest(config=config, atmosphere=atmosphere, star=star, windows=windows)


def test_durable_submission_reopens_saved_result(tmp_path, monkeypatch, workflow_daemon):
    from tumult.synthesis import operations, submission
    request = _request(Atmosphere('one.mod'), tmp_path)

    def synthesize(request, background_lines=True, cancel_token=None):
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), (None,), request, {}, {})

    monkeypatch.setattr(operations, 'synthesize', synthesize)
    handle = submission.submit_operation('synthesis.synthesize', request,
                                         root=tmp_path / 'jobs')
    result = handle.result(timeout=5)
    reopened = submission.open_workflow(handle.folder)
    assert reopened.status == 'finished'
    assert np.allclose(reopened.result().spectrum.values[1], result.spectrum.values[1])
    assert reopened.result().request.to_dict() == request.to_dict()
    assert reopened.progress['phase'] == 'finished'


def test_submission_records_scheduled_runs_before_completion(tmp_path, monkeypatch,
                                                            workflow_daemon):
    """The workflow record names its native runs as the executor schedules them."""
    import threading
    from tumult.synthesis import operations, runtime, submission
    request = _request(Atmosphere('one.mod'), tmp_path)
    scheduled = tmp_path / 'out' / 'run_a'
    release = threading.Event()
    seen = {}

    def synthesize(request, background_lines=True, cancel_token=None):
        runtime.note_run(str(scheduled))            # what _run_batch does per job
        seen['during'] = submission.open_workflow(handle.folder).runs
        release.wait(5)
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), (None,), request, {}, {}, runs=(str(scheduled),))

    monkeypatch.setattr(operations, 'synthesize', synthesize)
    handle = submission.submit_operation('synthesis.synthesize', request,
                                         root=tmp_path / 'jobs')
    release.set()
    handle.result(timeout=5)
    assert seen['during'] == (str(scheduled),)
    assert submission.open_workflow(handle.folder).runs == (str(scheduled),)
