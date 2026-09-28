"""The Slurm driver, and tool calculations recorded as workflows on the daemon."""
import json
import os
import time

import numpy as np
import pytest

from tumult.synthesis import Atmosphere, SynthConfig, SynthRequest, submission
from tumult.synthesis.cancel import CancelledError

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='detached workflows need POSIX')


def _request(tmp_path, **preferences):
    config = SynthConfig(root=tmp_path, cache='out', preferences=preferences)
    return SynthRequest(config=config, atmosphere=Atmosphere('one.mod'),
                        star=(5750, 4.5, 0, 1))


def _wait(predicate, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.1)
    return False


def test_slurm_driver_is_submitted_polled_and_cancelled(tmp_path, monkeypatch):
    from tumult.synthesis import scheduler
    calls = []
    live = {'state': 'queued'}
    monkeypatch.setattr(scheduler, 'sbatch', lambda script, args, env=None: (
        calls.append(('sbatch', str(script), list(args))) or '4242'))
    monkeypatch.setattr(scheduler, 'state', lambda job_id: live['state'])
    monkeypatch.setattr(scheduler, 'scancel', lambda job_id: (
        calls.append(('scancel', job_id)) or live.update(state=None) or True))
    script = tmp_path / 'driver.sbatch'
    script.write_text('#!/bin/sh\npython -m tumult.synthesis.driver "$1"\n')
    request = _request(tmp_path, workflow_runner='slurm', slurm_script='driver.sbatch')
    handle = submission.submit_operation('synthesis.synthesize', request,
                                         root=tmp_path / 'jobs')
    assert calls[0] == ('sbatch', str(script.resolve()), [str(handle.folder)])
    assert handle.to_dict()['driver']['slurm_id'] == '4242'
    assert handle.status == 'queued' and not handle.done
    live['state'] = 'running'
    assert handle.status == 'running'
    live['state'] = 'queued'
    assert handle.cancel()
    assert calls[-1] == ('scancel', '4242')
    assert handle.status == 'cancelled' and handle.done   # never started: settled here


def test_missing_slurm_script_fails_the_record_before_submission(tmp_path):
    with pytest.raises(ValueError):
        submission.submit_operation('synthesis.synthesize', _request(tmp_path),
                                    root=tmp_path / 'jobs', detach='slurm')
    records = list((tmp_path / 'jobs' / '_workflows').glob('*/workflow.json'))
    assert len(records) == 1 and json.loads(records[0].read_text())['state'] == 'failed'


def test_the_driver_runs_a_record_on_the_node_and_settles_it(tmp_path, monkeypatch):
    """What the sbatch script runs: the driver claims the record with its pid and the
    Slurm job id, executes it and exits 0 with the record finished."""
    from tumult import Spectrum
    from tumult.synthesis import driver, operations, tasks
    from tumult.synthesis.request import SynthesisResult

    def synthesize(request, background_lines=True, cancel_token=None):
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), (None,), request, {}, {})

    monkeypatch.setattr(operations, 'synthesize', synthesize)
    monkeypatch.setenv('SLURM_JOB_ID', '77')
    record = submission._new_record('node', None)
    record['task'] = {'kind': 'synthesis.synthesize',
                      'args': tasks.Operation.encode_args(_request(tmp_path), {})}
    record['driver'] = {'kind': 'slurm', 'slurm_id': None}
    folder = tmp_path / 'jobs' / '_workflows' / record['id']
    submission._atomic(folder / submission.RECORD, record)
    assert driver.main([str(folder)]) == 0
    settled = submission._read(folder)
    assert settled['state'] == 'finished' and settled['driver']['slurm_id'] == '77'
    assert settled['driver']['pid'] == os.getpid()


# --- task records: the tool calculations on the daemon -------------------------------------
# A fit record carries the request and the fit's keyword arguments; its stub sees the
# decoded objects and returns a SpectrumFit.
_FAKE_FIT = '''
import time
from pathlib import Path
import numpy as np
from tumult.synthesis import operations
from tumult.synthesis.request import SpectrumFit
from tumult import Spectrum

def fit(request, *, target, cancel_token=None, **kwargs):
    release = Path(request.config.preferences['release_file'])
    while not release.exists():
        cancel_token.check()
        time.sleep(0.05)
    x, y = target.values
    return SpectrumFit({'teff': 5600.0, 'vshift': 1.5}, Spectrum(x, y * 0.5), {'vshift': 1.5},
                       request)

operations.fit = fit
'''


@pytest.fixture
def fit_env(tmp_path, monkeypatch, daemon_dir):
    """Spawned daemons run the fake fit in this test's own daemon directory."""
    hook = tmp_path / 'hook'
    hook.mkdir()
    (hook / 'sitecustomize.py').write_text(_FAKE_FIT)
    monkeypatch.setenv('PYTHONPATH', str(hook))
    return tmp_path / 'release'


def test_operation_records_run_on_the_daemon_and_load_their_kind(tmp_path, fit_env):
    """A fit submitted as an operation record is claimed by the daemon, decoded into the
    live objects the calculation takes, and its result reopens as a SpectrumFit from any
    process; the record's task block names the operation."""
    from tumult import Spectrum
    request = _request(tmp_path, release_file=str(fit_env)).update(windows=[(5000., 5001.)])
    obs = Spectrum(np.array([5000., 5000.5, 5001.]), np.array([1., .8, 1.]))
    handle = submission.submit_operation('synthesis.fit', request, root=tmp_path / 'jobs',
                                         label='fit obs', target=obs,
                                         fit_post={'vshift': 0.0})
    assert handle.kind.name == 'synthesis.fit' and handle.label == 'fit obs'
    record = handle.to_dict()
    assert record['task']['kind'] == 'synthesis.fit'
    assert record['task']['args']['kwargs']['fit_post'] == {'vshift': 0.0}
    assert _wait(lambda: handle.status == 'running')
    fit_env.write_text('')
    result = submission.open_workflow(handle.folder).result(timeout=20)
    assert result.best == {'teff': 5600.0, 'vshift': 1.5} and result.post == {'vshift': 1.5}
    assert np.allclose(result.spectrum.values[1], [.5, .4, .5])
    assert result.request.to_dict() == request.to_dict()
    assert handle.status == 'finished'


def test_task_progress_snapshots_are_written_and_read_back(tmp_path, monkeypatch,
                                                           workflow_daemon):
    """A streaming kind's intermediate results reach the record's progress file, at
    most every PROGRESS_MIN_S, and read back decoded; a kind without progress has none.
    Runs on a daemon in this process with a stubbed balance solver."""
    from tumult import sparams
    from tumult.measurement import BalanceResult

    def ionbalance(request, *, measured, progress=None, cancel_token=None, **kwargs):
        result = measured
        for n in range(3):
            result.history.append({'params': {'teff': 5000. + n}})
            progress(result)
        result.stage = 'solved'
        return result

    monkeypatch.setattr(sparams, 'ionbalance', ionbalance)
    monkeypatch.setattr(submission, 'PROGRESS_MIN_S', 0.0)
    from tumult.synthesis import daemon
    events = []
    listener = daemon.Client('client')
    listener.listen(events.append)
    result = BalanceResult(element='Fe', uid='bal-1', params0={'teff': 5000.})
    handle = submission.submit_operation('sparams.ionbalance', _request(tmp_path),
                                         measured=result, max_iter=3)
    assert handle.folder.parent == (tmp_path / 'out' / '_workflows').resolve()
    solved = handle.result(timeout=10)
    assert isinstance(solved, BalanceResult) and solved.uid == 'bal-1'
    assert solved.stage == 'solved' and len(solved.history) == 3
    snapshot = handle.progress_result()
    assert isinstance(snapshot, BalanceResult) and handle.progress_stamp() is not None
    assert 1 <= len(snapshot.history) <= 3            # whichever round was written last
    assert _wait(lambda: [e['stamp'] for e in events if e['event'] == 'progress'
                          and e['job'] == handle.id][-1:] == [handle.progress_stamp()], 5)
    listener.close()

    with pytest.raises(KeyError, match='unknown workflow kind'):
        from tumult.synthesis import tasks
        tasks.kind_of({'task': {'kind': 'no-such-kind'}})
