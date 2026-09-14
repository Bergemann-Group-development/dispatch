"""The workflow daemon: lifetime rules, reconnection, latency and thread accounting.

Lifetime tests run the daemon as the process a client spawns. Its calculation is
replaced there through a ``sitecustomize`` on its ``PYTHONPATH``: a synthesis waits for
the request's ``release_file`` or sleeps ``sleep`` seconds. Tests that need to see inside
the calculation run the daemon in a thread instead (the ``workflow_daemon`` fixture).
"""
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from tumult.synthesis import Atmosphere1D, NativeRequest, daemon, submission
from tumult.synthesis.cancel import CancelledError

pytestmark = pytest.mark.skipif(os.name != 'posix', reason='the daemon needs POSIX')

ROOT = Path(__file__).resolve().parent.parent

_FAKE_EXECUTE = '''
import time
from pathlib import Path
import numpy as np
from tumult.synthesis import workflows
from tumult.synthesis.request import SynthesisResult
from tumult import Spectrum

def execute(resolved, cancel_token=None):
    prefs = resolved.request.prefs
    if prefs.get('sleep'):
        time.sleep(float(prefs['sleep']))
    release = Path(prefs.get('release_file') or '/nonexistent')
    while prefs.get('release_file') and not release.exists():
        cancel_token.check()
        time.sleep(0.02)
    spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
    return SynthesisResult((spec,), ({'window': None},), resolved.request, resolved.settings)

workflows.execute_request = execute
'''


@pytest.fixture
def fake_env(tmp_path, monkeypatch):
    """Spawned daemons run the fake calculation; returns the release file."""
    hook = tmp_path / 'hook'
    hook.mkdir()
    (hook / 'sitecustomize.py').write_text(_FAKE_EXECUTE)
    monkeypatch.setenv('PYTHONPATH', str(hook))
    return tmp_path / 'release'


def _request(tmp_path, **prefs):
    return NativeRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out', **prefs},
                        run_vals={'mode': 'LTE'},
                        selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))


def _wait(predicate, timeout=20, step=0.02):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(step)
    return False


def _daemon_pid():
    info = daemon.running()
    return info.get('pid') if info else None


_SERVER = '''
import sys, time
sys.path.insert(0, {root!r})
from tumult.synthesis import daemon, submission, NativeRequest, Atmosphere1D
server = daemon.connect_server()
request = NativeRequest({prefs!r}, run_vals={{'mode': 'LTE'}},
                       selection=Atmosphere1D('one.mod', (5750, 4.5, 0, 1)))
handle = submission.submit(request, root={jobs!r}) if {submit!r} else None
print(handle.id if handle else 'connected', flush=True)
time.sleep(600)
'''


def _server_process(tmp_path, submit, **prefs):
    """A server process connected to the daemon, optionally having submitted a job;
    returns the process and the line it printed (the workflow id)."""
    script = _SERVER.format(root=str(ROOT), jobs=str(tmp_path / 'jobs'), submit=submit,
                            prefs={'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out',
                                   **prefs})
    proc = subprocess.Popen([sys.executable, '-c', script], stdout=subprocess.PIPE, text=True)
    line = proc.stdout.readline().strip()
    return proc, line


def test_rule_a_a_killed_server_leaves_its_job_to_finish_then_the_daemon_exits(tmp_path,
                                                                               fake_env):
    proc, workflow_id = _server_process(tmp_path, True, release_file=str(fake_env))
    handle = submission.open_workflow(tmp_path / 'jobs' / '_workflows' / workflow_id)
    assert _wait(lambda: handle.status == 'running')
    pid = _daemon_pid()
    proc.kill()
    proc.wait()
    time.sleep(1.0)
    assert _daemon_pid() == pid                     # still running the job
    fake_env.write_text('')
    assert _wait(lambda: handle.status == 'finished')
    assert _wait(lambda: daemon.running() is None, timeout=10)


def test_rule_b_a_server_disconnecting_from_an_idle_daemon_ends_it(tmp_path, fake_env):
    server = daemon.connect_server()
    assert _daemon_pid() == server.info['daemon_pid']
    time.sleep(daemon.RECENT_S)
    daemon.disconnect()
    started = time.monotonic()
    assert _wait(lambda: daemon.running() is None, timeout=5)
    assert time.monotonic() - started < 1.0


def test_a_client_does_not_keep_an_idle_daemon_alive_but_its_work_does(tmp_path, fake_env):
    client = daemon.client()
    assert client.info['daemon_pid']
    assert _wait(lambda: daemon.running() is None, timeout=5)     # connection still open
    assert not client.connected or _wait(lambda: not client.connected, timeout=2)

    handle = submission.submit(_request(tmp_path, release_file=str(fake_env)),
                               root=tmp_path / 'jobs')
    daemon.disconnect()                                           # the script exits
    assert _wait(lambda: handle.status == 'running')
    time.sleep(1.5)
    assert daemon.running() is not None
    fake_env.write_text('')
    assert _wait(lambda: handle.status == 'finished')
    assert _wait(lambda: daemon.running() is None, timeout=10)


def test_two_servers_starting_at_once_share_one_daemon(tmp_path, fake_env):
    barrier = threading.Barrier(2)
    clients = []

    def start():
        barrier.wait()
        clients.append(daemon.Client('server'))

    threads = [threading.Thread(target=start) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(30)
    try:
        assert len(clients) == 2
        assert clients[0].info['daemon_pid'] == clients[1].info['daemon_pid'] == _daemon_pid()
    finally:
        for client in clients:
            client.close()


def test_a_reconnecting_server_lists_a_running_job_and_receives_its_settle_event(tmp_path,
                                                                                fake_env):
    first = daemon.connect_server()
    handle = submission.submit(_request(tmp_path, release_file=str(fake_env)),
                               root=tmp_path / 'jobs')
    assert _wait(lambda: handle.status == 'running')
    pid = first.info['daemon_pid']
    daemon.disconnect()

    events = []
    second = daemon.connect_server(roots=[tmp_path / 'out'])
    assert second.info['daemon_pid'] == pid
    listed = {job['job']: job for job in second.request('list')['jobs']}
    assert listed[handle.id]['state'] == 'running' and listed[handle.id]['label'] is None
    second.follow(handle.id, events.append)
    fake_env.write_text('')
    assert _wait(lambda: any(e.get('state') == 'finished' for e in events))
    assert handle.status == 'finished'


def test_a_short_job_returns_quickly_on_a_warm_daemon(tmp_path, fake_env):
    daemon.connect_server()
    submission.submit(_request(tmp_path, sleep=0.1), root=tmp_path / 'jobs').result(timeout=30)
    request = _request(tmp_path, sleep=0.1)
    started = time.monotonic()
    handle = submission.submit(request, root=tmp_path / 'jobs')
    assert np.allclose(handle.result(timeout=10).spectrum.values[1], [1., .9])
    assert time.monotonic() - started < 0.3


def test_a_daemon_killed_mid_job_leaves_an_interrupted_record_for_the_next(tmp_path,
                                                                          fake_env):
    daemon.client()
    handle = submission.submit(_request(tmp_path, release_file=str(fake_env)),
                               root=tmp_path / 'jobs')
    assert _wait(lambda: handle.status == 'running')
    pid = _daemon_pid()
    os.kill(pid, signal.SIGKILL)
    assert _wait(lambda: handle.status == 'interrupted')          # no daemon holds it
    assert handle.to_dict()['state'] == 'interrupted'
    daemon.disconnect()
    next_daemon = daemon.client()                                 # registers roots.json
    assert next_daemon.info['daemon_pid'] != pid
    record = json.loads((handle.folder / submission.RECORD).read_text())
    assert record['state'] == 'interrupted' and record['error']
    with pytest.raises(RuntimeError, match='interrupted'):
        handle.result()


def test_a_queued_record_of_a_dead_daemon_is_run_by_the_next(tmp_path, fake_env):
    folder = tmp_path / 'jobs' / '_workflows' / 'w1'
    record = submission._new_record('orphan', None)
    record.update(id='w1', driver={'kind': 'daemon', 'host': None, 'pid': None},
                  task={'kind': 'synthesis'})
    from tumult.synthesis.workflows import resolve_request
    record['request'] = resolve_request(_request(tmp_path, sleep=0.05)).to_dict()
    submission._atomic(folder / submission.RECORD, record)
    handle = submission.open_workflow(folder)
    assert handle.status == 'queued' and not handle.done
    assert handle.result(timeout=30) is not None           # the handle hands it over


# --- in a thread: what the daemon runs is visible to the test ---------------------------
def _stub(monkeypatch, execute):
    from tumult.synthesis import workflows
    monkeypatch.setattr(workflows, 'execute_request', execute)


def _result(resolved):
    from tumult import Spectrum
    from tumult.synthesis.request import SynthesisResult
    spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
    return SynthesisResult((spec,), ({'window': None},), resolved.request, resolved.settings)


def test_workflows_run_one_at_a_time_in_submission_order(tmp_path, monkeypatch,
                                                         workflow_daemon):
    order, release = [], threading.Event()

    def execute(resolved, cancel_token=None):
        order.append(resolved.request.prefs['n'])
        assert release.wait(10)
        return _result(resolved)

    _stub(monkeypatch, execute)
    handles = [submission.submit(_request(tmp_path, n=n), root=tmp_path / 'jobs')
               for n in range(3)]
    assert _wait(lambda: handles[0].status == 'running')
    assert [h.status for h in handles[1:]] == ['queued', 'queued']
    release.set()
    for handle in handles:
        handle.result(timeout=10)
    assert order == [0, 1, 2]


def test_cancelling_queued_and_running_workflows(tmp_path, monkeypatch, workflow_daemon):
    started = threading.Event()

    def execute(resolved, cancel_token=None):
        started.set()
        while True:
            cancel_token.check()
            time.sleep(0.01)

    _stub(monkeypatch, execute)
    running = submission.submit(_request(tmp_path), root=tmp_path / 'jobs')
    queued = submission.submit(_request(tmp_path), root=tmp_path / 'jobs')
    assert started.wait(10)
    seen = []
    queued.add_done_callback(lambda h: seen.append(h.status))
    assert submission.open_workflow(queued.folder).cancel()
    assert queued.status == 'cancelled' and not queued.to_dict()['driver']['pid']
    assert _wait(lambda: seen == ['cancelled'], timeout=2)
    assert running.cancel()
    with pytest.raises(CancelledError):
        running.result(timeout=10)
    assert not running.cancel()                 # nothing left to cancel


def test_thread_events_follow_the_allocation_and_a_shared_run_counts_once(
        tmp_path, monkeypatch, workflow_daemon):
    """Two batches of one workflow ask for the same output: the service runs it once, the
    workflow's ``threads`` event reports that run's allocation, and ``usage`` matches."""
    import contextvars
    from concurrent.futures import ThreadPoolExecutor
    from tumult import launch
    from tumult.synthesis.executor import _run_batch
    from tumult.synthesis.runtime import get_runtime
    import tumult as tm

    executable = tmp_path / 'dispatch.x'
    gate = tmp_path / 'gate'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys, time\nfrom pathlib import Path\n'
        'p=Path(sys.argv[1]).with_suffix(""); p.mkdir(parents=True, exist_ok=True)\n'
        f'while not Path({str(gate)!r}).exists(): time.sleep(.01)\n'
        '(p / "abund.txt").write_text("finished")\n')
    executable.chmod(0o755)
    monkeypatch.setattr(tm, 'read', lambda path: path.name)
    budget = get_runtime().budget
    allocations = []

    def execute(resolved, cancel_token=None):
        jobs = [('a', tmp_path / 'shared', {})]
        kwargs = {'cwd': tmp_path, 'verbose': 0}
        with ThreadPoolExecutor(2) as pool:
            futures = [pool.submit(contextvars.copy_context().run, _run_batch, jobs, kwargs,
                                   budget, cancel_token) for _ in range(2)]
            assert _wait(lambda: get_runtime().used > 0, timeout=10)
            _wait(lambda: False, timeout=0.3)
            with get_runtime()._lock:
                allocations.append(sum(t.allocated for t in get_runtime()._running))
            gate.write_text('')
            for future in futures:
                future.result(10)
        return _result(resolved)

    _stub(monkeypatch, execute)
    events = []
    listener = daemon.Client('client')
    listener.listen(events.append)
    handle = submission.submit(_request(tmp_path), root=tmp_path / 'jobs')
    handle.result(timeout=30)
    assert _wait(lambda: any(e['event'] == 'threads' and e['threads'] == 0 for e in events))
    listener.close()
    threads = [e['threads'] for e in events
               if e['event'] == 'threads' and e['job'] == handle.id]
    usage = [e['used'] for e in events if e['event'] == 'usage']
    assert allocations == [budget]              # one run holds the whole budget, not two
    assert max(threads) == budget and threads[-1] == 0
    assert max(usage) == budget and usage[-1] == 0


def test_the_list_names_live_external_runs_and_cancels_them(tmp_path, monkeypatch,
                                                            workflow_daemon):
    import socket
    from tumult.synthesis import jobrecord
    root = tmp_path / 'out'
    run = root / 'shell_run'
    run.mkdir(parents=True)
    (run / 'job.nml').write_text(
        "&job_record\n state='running', pid=%d, host='%s', input='shell_run.nml', "
        "omp_threads=4, mpi_size=2,\n /\n" % (os.getpid(), socket.gethostname()))
    finished = root / 'done_run'
    finished.mkdir()
    (finished / 'job.nml').write_text("&job_record\n state='finished', pid=1,\n /\n")
    server = daemon.connect_server(roots=[root])
    workflow_daemon.scan_external()
    external = [job for job in server.request('list')['jobs'] if job['type'] == 'external']
    assert [(job['label'], job['threads']) for job in external] == [('shell_run', 8)]
    signalled = []
    monkeypatch.setattr(jobrecord, 'cancel', lambda folder: signalled.append(folder) or True)
    assert server.request('cancel', job=external[0]['job'])['state'] == 'cancelled'
    assert signalled == [str(run)]
