"""The Slurm driver, and tool calculations recorded as workflows on the daemon."""
import json
import os
import time

import numpy as np
import pytest

from tumult.synthesis import Atmosphere1D, NativeRequest, submission
from tumult.synthesis.cancel import CancelledError

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='detached workflows need POSIX')


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
    request = NativeRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
                            'workflow_runner': 'slurm', 'slurm_script': 'driver.sbatch'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    handle = submission.submit(request, root=tmp_path / 'jobs')
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
    request = NativeRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    with pytest.raises(ValueError):
        submission.submit(request, root=tmp_path / 'jobs', detach='slurm')
    records = list((tmp_path / 'jobs' / '_workflows').glob('*/workflow.json'))
    assert len(records) == 1 and json.loads(records[0].read_text())['state'] == 'failed'


def test_the_driver_runs_a_record_on_the_node_and_settles_it(tmp_path, monkeypatch):
    """What the sbatch script runs: the driver claims the record with its pid and the
    Slurm job id, executes it and exits 0 with the record finished."""
    from tumult import Spectrum
    from tumult.synthesis import driver, workflows
    from tumult.synthesis.request import SynthesisResult

    def execute(resolved, cancel_token=None):
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), ({'window': None},), resolved.request, resolved.settings)

    monkeypatch.setattr(workflows, 'execute_request', execute)
    monkeypatch.setenv('SLURM_JOB_ID', '77')
    request = NativeRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    record = submission._new_record('node', None)
    record['request'] = workflows.resolve_request(request).to_dict()
    record['driver'] = {'kind': 'slurm', 'slurm_id': None}
    folder = tmp_path / 'jobs' / '_workflows' / record['id']
    submission._atomic(folder / submission.RECORD, record)
    assert driver.main([str(folder)]) == 0
    settled = submission._read(folder)
    assert settled['state'] == 'finished' and settled['driver']['slurm_id'] == '77'
    assert settled['driver']['pid'] == os.getpid()


# --- task records: the tool calculations on the daemon -------------------------------------
# A fit record carries the observation and the FitRequest; its stub sees the decoded
# objects, reports one progress snapshot through the fit kind's result encoding when the
# kind streams, and returns (best, raw, post).
_FAKE_FIT = '''
import time
from pathlib import Path
import numpy as np
from tumult.synthesis import workflows, tasks
from tumult import Spectrum

def fit(observation, request, *, cancel_token=None):
    release = Path(request.synthesis.prefs['release_file'])
    while not release.exists():
        cancel_token.check()
        time.sleep(0.05)
    x, y = observation.values
    return ({'teff': 5600.0, 'vshift': 1.5}, Spectrum(x, y * 0.5), {'vshift': 1.5})

workflows.fit = fit
'''


@pytest.fixture
def fit_env(tmp_path, monkeypatch):
    hook = tmp_path / 'hook'
    hook.mkdir()
    (hook / 'sitecustomize.py').write_text(_FAKE_FIT)
    monkeypatch.setenv('PYTHONPATH', str(hook))
    return tmp_path / 'release'


def _fit_args(tmp_path, release):
    from tumult import Spectrum
    from tumult.synthesis import FitRequest, MaskWindow
    prefs = {'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
             'release_file': str(release)}
    request = NativeRequest(prefs, run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)),
                           masks=(MaskWindow((5000., 5001.), (5000., 5001.)),))
    fit_request = FitRequest(request, fit_post={'vshift': 0.0})
    obs = Spectrum(np.array([5000., 5000.5, 5001.]), np.array([1., .8, 1.]))
    return {'observation': obs.to_dict(), 'request': fit_request.to_dict()}, prefs


def test_task_records_run_on_the_daemon_and_load_their_kind(tmp_path, fit_env):
    """A fit submitted as a task record is claimed by the daemon, decoded into the live
    objects the calculation takes, and its result reopens as (best, raw, post) from any
    process; the record's task block names the kind."""
    args, prefs = _fit_args(tmp_path, fit_env)
    handle = submission.submit_task('fit', args, prefs=prefs, root=tmp_path / 'jobs',
                                    label='fit obs')
    assert handle.kind.name == 'fit' and handle.label == 'fit obs'
    record = handle.to_dict()
    assert record['task']['kind'] == 'fit' and record['request'] is None
    assert record['task']['args']['request']['fit_post'] == {'vshift': 0.0}
    assert _wait(lambda: handle.status == 'running')
    fit_env.write_text('')
    best, raw, post = submission.open_workflow(handle.folder).result(timeout=20)
    assert best == {'teff': 5600.0, 'vshift': 1.5} and post == {'vshift': 1.5}
    assert np.allclose(raw.values[1], [.5, .4, .5])
    assert handle.status == 'finished'


def test_task_progress_snapshots_are_written_and_read_back(tmp_path, monkeypatch,
                                                           workflow_daemon):
    """A streaming kind's intermediate results reach the record's progress file, at
    most every PROGRESS_MIN_S, and read back decoded; a kind without progress has none.
    Runs on a daemon in this process with a stubbed balance solver."""
    from tumult import sparams
    from tumult.measurement import BalanceResult

    def solve(request, result, opts=None, *, progress=None, cancel_token=None):
        for n in range(3):
            result.history.append({'params': {'teff': 5000. + n}})
            progress(result)
        result.stage = 'solved'
        return result

    monkeypatch.setattr(sparams, 'solve_ionisation_balance', solve)
    monkeypatch.setattr(submission, 'PROGRESS_MIN_S', 0.0)
    from tumult.synthesis import daemon
    events = []
    listener = daemon.Client('client')
    listener.listen(events.append)
    prefs = {'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'}
    request = NativeRequest(prefs, run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    result = BalanceResult(element='Fe', uid='bal-1', params0={'teff': 5000.})
    args = {'request': request.to_dict(), 'result': result.to_dict(), 'masks': None,
            'observation': None, 'linelist': None, 'opts': {'max_iter': 3}}
    handle = submission.submit_task('ionisation_balance', args, prefs=prefs)
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
