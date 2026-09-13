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


# --- the daemon --------------------------------------------------------------------------
@pytest.fixture
def stop_daemon(tmp_path):
    """Whatever daemon a test started for tmp_path is ended afterwards."""
    from tumult.synthesis import daemon
    yield
    live = daemon.status(tmp_path / 'jobs')
    if live is not None:
        os.kill(live[0], 15)
        _wait(lambda: daemon.status(tmp_path / 'jobs') is None)


def test_daemon_runs_workflows_off_the_server_and_is_reused(tmp_path, driver_env, stop_daemon):
    from tumult.synthesis import daemon
    first = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                              detach='daemon')
    assert first.detached == 'daemon'
    assert _wait(lambda: daemon.status(tmp_path / 'jobs') is not None)
    pid, _host = daemon.status(tmp_path / 'jobs')
    assert _wait(lambda: first.status == 'running')
    assert first.to_dict()['driver']['pid'] == pid
    second = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                               detach='daemon')
    assert _wait(lambda: second.status == 'running')
    assert daemon.status(tmp_path / 'jobs')[0] == pid          # one daemon, two workflows
    driver_env.write_text('')
    assert np.allclose(submission.open_workflow(first.folder).result(timeout=20)
                       .spectrum.values[1], [1., .9])
    assert second.result(timeout=20) is not None
    assert first.status == second.status == 'finished'


def test_daemon_honours_a_cancel_flag_from_any_handle(tmp_path, driver_env, stop_daemon):
    handle = submission.submit(_request(tmp_path, driver_env), root=tmp_path / 'jobs',
                               detach='daemon')
    assert _wait(lambda: handle.status == 'running')
    seen = []
    handle.add_done_callback(lambda h: seen.append(h.status))
    assert submission.open_workflow(handle.folder).cancel()
    assert handle.to_dict().get('cancel_requested') is True
    assert _wait(lambda: handle.status == 'cancelled')
    with pytest.raises(CancelledError):
        handle.result(timeout=5)
    assert _wait(lambda: seen == ['cancelled'], timeout=3 * submission.POLL_S)


def test_daemon_records_wait_queued_without_a_daemon_and_are_revived(tmp_path, monkeypatch):
    """With no daemon, a queued record stays queued (the next daemon claims it), a claimed
    one is interrupted, and cancelling a queued one settles it in place."""
    from tumult.synthesis import daemon
    monkeypatch.setattr(daemon, 'status', lambda root: None)
    started = []
    monkeypatch.setattr(daemon, 'ensure', lambda root: started.append(root))
    request = SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    handle = submission.submit(request, root=tmp_path / 'jobs', detach='daemon')
    assert started == [str((tmp_path / 'jobs').resolve())]
    assert handle.status == 'queued' and not handle.done
    # a reconnecting server finds the queued work and starts a daemon for it
    submission.discover_workflows(tmp_path / 'jobs')
    assert len(started) == 2
    submission.mutate_record(handle.folder, lambda r: r['driver'].update(pid=4321))
    submission.update_record(handle.folder, state='running')
    assert handle.status == 'interrupted' and not handle.cancel()
    submission.mutate_record(handle.folder, lambda r: r['driver'].update(pid=None))
    submission.update_record(handle.folder, state='queued')
    assert handle.cancel() and handle.status == 'cancelled'


def test_daemon_exits_when_idle_and_claims_once(tmp_path, monkeypatch):
    """The daemon class in a thread: it claims a queued record exactly once, runs it on
    the WORKFLOWS pool with the stubbed execution, and leaves after the idle period."""
    import threading
    from tumult.synthesis import daemon, workflows
    from tumult.synthesis.request import SynthesisResult
    from tumult import Spectrum

    def execute(resolved, cancel_token=None):
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), ({'window': None},), resolved.request, resolved.settings)

    monkeypatch.setattr(workflows, 'execute_request', execute)
    monkeypatch.setattr(daemon, 'ensure', lambda root: None)
    request = SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out'},
                           run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
    handle = submission.submit(request, root=tmp_path / 'jobs', detach='daemon')
    worker = daemon.Daemon(tmp_path / 'jobs', idle_s=1.0)
    thread = threading.Thread(target=worker.serve, daemon=True)
    thread.start()
    assert _wait(lambda: (tmp_path / 'jobs' / '_workflows' / daemon.PIDFILE).is_file())
    assert np.allclose(handle.result(timeout=10).spectrum.values[1], [1., .9])
    assert handle.to_dict()['driver']['pid'] == os.getpid()
    thread.join(10)
    assert not thread.is_alive()
    assert not (tmp_path / 'jobs' / '_workflows' / daemon.PIDFILE).exists()


# --- task records: the tool calculations on the same drivers ------------------------------
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
    request = SynthRequest(prefs, run_vals={'mode': 'LTE'},
                           selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)),
                           masks=(MaskWindow((5000., 5001.), (5000., 5001.)),))
    fit_request = FitRequest(request, fit_post={'vshift': 0.0})
    obs = Spectrum(np.array([5000., 5000.5, 5001.]), np.array([1., .8, 1.]))
    return {'observation': obs.to_dict(), 'request': fit_request.to_dict()}, prefs


def test_task_records_run_on_the_daemon_and_load_their_kind(tmp_path, fit_env, stop_daemon):
    """A fit submitted as a task record is claimed by the daemon, decoded into the live
    objects the calculation takes, and its result reopens as (best, raw, post) from any
    process; the record's task block names the kind."""
    args, prefs = _fit_args(tmp_path, fit_env)
    handle = submission.submit_task('fit', args, prefs=prefs, root=tmp_path / 'jobs',
                                    label='fit obs', detach='daemon')
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


def test_task_progress_snapshots_are_written_and_read_back(tmp_path, monkeypatch):
    """A streaming kind's intermediate results reach the record's progress file, at
    most every PROGRESS_MIN_S, and read back decoded; a kind without progress has none.
    Runs in this process with a stubbed balance solver."""
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
    prefs = {'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
             'workflow_runner': 'in-process'}
    request = SynthRequest(prefs, run_vals={'mode': 'LTE'},
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

    with pytest.raises(KeyError, match='unknown workflow kind'):
        from tumult.synthesis import tasks
        tasks.kind_of({'task': {'kind': 'no-such-kind'}})
