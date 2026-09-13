"""Detached workflow drivers: a process that outlives the submitter, and Slurm."""
import json
import os
import time

import numpy as np
import pytest

from tumult.synthesis import Atmosphere1D, SynthRequest, submission
from tumult.synthesis.cancel import CancelledError

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='detached drivers need POSIX')

# The driver runs in a fresh interpreter, so the calculation is replaced there through a
# sitecustomize on its PYTHONPATH rather than by monkeypatching this process.
_FAKE_EXECUTE = '''
import time
from pathlib import Path
import numpy as np
from tumult.synthesis import workflows
from tumult.synthesis.request import SynthesisResult
from tumult import Spectrum

def execute(resolved, cancel_token=None):
    release = Path(resolved.request.prefs['release_file'])
    while not release.exists():
        cancel_token.check()
        time.sleep(0.05)
    spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
    return SynthesisResult((spec,), ({'window': None},), resolved.request, resolved.settings)

workflows.execute_request = execute
'''


@pytest.fixture
def driver_env(tmp_path, monkeypatch):
    hook = tmp_path / 'hook'
    hook.mkdir()
    (hook / 'sitecustomize.py').write_text(_FAKE_EXECUTE)
    monkeypatch.setenv('PYTHONPATH', str(hook))
    return tmp_path / 'release'


def _request(tmp_path, release):
    return SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
                         'release_file': str(release)},
                        run_vals={'mode': 'LTE'},
                        selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))


def _wait(predicate, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.1)
    return False


def test_process_driver_outlives_the_submitter_and_delivers(tmp_path, driver_env):
    handle = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                               detach='process')
    assert handle.detached == 'process'
    assert _wait(lambda: handle.status == 'running')
    record = handle.to_dict()
    assert record['driver']['pid'] and record['driver']['host']
    assert not handle.done
    # a second process's view: nothing of the submitter is needed
    reopened = submission.open_workflow(handle.folder)
    assert reopened.status == 'running' and not reopened.done
    driver_env.write_text('')
    result = reopened.result(timeout=20)
    assert np.allclose(result.spectrum.values[1], [1., .9])
    assert handle.status == 'finished'
    assert str(handle.folder / submission.DRIVER_LOG) in handle.logs


def test_process_driver_is_cancelled_by_signal_and_settles_its_record(tmp_path, driver_env):
    handle = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                               detach='process')
    # wait for the driver to claim the record, so the signal reaches its handler and the
    # cancellation runs through the driver's own token rather than a startup kill
    assert _wait(lambda: (handle.to_dict()['driver'] or {}).get('started'))
    seen = []
    handle.add_done_callback(lambda h: seen.append(h.status))
    assert submission.open_workflow(handle.folder).cancel()
    assert _wait(lambda: handle.status == 'cancelled')
    with pytest.raises(CancelledError):
        handle.result(timeout=5)
    assert _wait(lambda: seen == ['cancelled'], timeout=3 * submission.POLL_S)
    assert not handle.cancel()          # nothing left to signal
    assert 'cancelling' in (handle.folder / submission.DRIVER_LOG).read_text()


def test_process_driver_killed_during_startup_is_settled_by_the_canceller(tmp_path, driver_env):
    handle = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                               detach='process')
    assert _wait(lambda: handle.status == 'running')      # the pid is alive, unclaimed
    assert handle.cancel()
    assert handle.status == 'cancelled'


def test_dead_driver_reports_interrupted_and_is_never_signalled(tmp_path):
    folder = tmp_path / 'jobs' / '_workflows' / 'w1'
    folder.mkdir(parents=True)
    # a pid that cannot be ours and whose command line cannot name this workflow
    (folder / submission.RECORD).write_text(json.dumps({
        'id': 'w1', 'state': 'running', 'phase': 'synthesis', 'submitted': time.time(),
        'runs': [], 'driver': {'kind': 'process', 'pid': os.getpid(), 'host': None}}))
    handle = submission.open_workflow(folder)
    assert handle.status == 'interrupted' and handle.done
    assert not handle.cancel()
    with pytest.raises(RuntimeError):
        handle.result()


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
    request = SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
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
    request = SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    with pytest.raises(ValueError):
        submission.submit(request, root=tmp_path / 'jobs', detach='slurm')
    records = list((tmp_path / 'jobs' / '_workflows').glob('*/workflow.json'))
    assert len(records) == 1 and json.loads(records[0].read_text())['state'] == 'failed'
