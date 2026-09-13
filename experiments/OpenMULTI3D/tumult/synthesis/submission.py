"""Durable asynchronous submission of the same resolved calculation used synchronously.

A workflow is recorded in ``<root>/_workflows/<id>/workflow.json`` and executed by one
of three drivers, chosen at submission:

* in-process (the default): a thread of this process's ``WORKFLOWS`` pool. It ends with
  the process; a server shutdown cancels it.
* ``daemon``: the long-lived worker of :mod:`tumult.synthesis.daemon` serving this root,
  started on demand. It outlives the submitter without a per-workflow interpreter start;
  a cancel is a ``cancel_requested`` flag in the record, which the daemon honours.
* ``process``: ``python -m tumult.synthesis.driver <folder>`` in its own session on this
  host. It outlives the submitter, and a later process cancels it by signalling the
  recorded process group once the pid is verified to still run that driver.
* ``slurm``: the same driver submitted through ``sbatch <script> <folder>``; the
  scheduler owns liveness (``squeue``) and cancellation (``scancel``).

The detached drivers settle the record themselves, so the record is shared between
processes and every update takes the folder's file lock. A handle opened on a detached
record reports its state from the driver's liveness; callbacks on such a handle are
delivered by a polling watcher rather than by a future.
"""
from __future__ import annotations

from copy import deepcopy
from contextlib import contextmanager
import json
import logging
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import threading
import time
import uuid

from .cancel import CancelToken, CancelledError
from .jobrecord import process_runs as _process_alive
from .request import ResolvedSynthRequest, SynthesisResult, SynthRequest
from .runtime import WORKFLOWS, observe_runs

RECORD = 'workflow.json'
LOCK = 'workflow.lock'
DRIVER_LOG = 'driver.log'
DETACHED = ('daemon', 'process', 'slurm')
# preference spellings of the driver choice (TumultStudio2 Preferences -> Engine)
RUNNER_PREFS = {'local': 'daemon', 'daemon': 'daemon', 'process': 'process',
                'slurm': 'slurm', 'in-process': None, '': None}
TERMINAL = ('finished', 'failed', 'cancelled', 'interrupted')
POLL_S = 3.0        # how often the watcher re-reads detached records for callbacks
_live = {}
_lock = threading.RLock()
log = logging.getLogger(__name__)


# --- the record on disk ----------------------------------------------------------------
def _atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'{path.name}.{os.getpid()}.partial')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + '\n')
    tmp.replace(path)


def _read(folder):
    try:
        return json.loads((Path(folder) / RECORD).read_text())
    except (OSError, ValueError):
        return None


@contextmanager
def _locked(folder):
    """Hold the record's advisory lock across a read-modify-write. The lock file lives
    beside the record so a driver in another process (or on another node sharing the
    file system) serialises with the submitter and with a cancelling server."""
    with _lock:
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        with open(folder / LOCK, 'a+') as handle:
            if os.name == 'posix':
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                if os.name == 'posix':
                    fcntl.flock(handle, fcntl.LOCK_UN)


def update_record(folder, base=None, **values):
    """Merge ``values`` into the record under its lock; ``base`` stands in for a record
    that cannot be read back. Returns the record as written."""
    with _locked(folder):
        current = _read(folder) or dict(base or {})
        current.update(values)
        _atomic(Path(folder) / RECORD, current)
        return current


def mutate_record(folder, change):
    """Apply ``change(record)`` in place under the lock, for nested fields."""
    with _locked(folder):
        current = _read(folder)
        if current is None:
            raise FileNotFoundError(f'No workflow record in {folder}')
        change(current)
        _atomic(Path(folder) / RECORD, current)
        return current


# --- the driver's liveness -------------------------------------------------------------
def driver_state(record):
    """The live state of a detached record's driver: ``queued``, ``running``, ``False``
    (gone) or ``None`` (not observable from here: an in-process record, a process driver
    on another host, a platform without ``ps``)."""
    driver = record.get('driver') or {}
    kind = driver.get('kind')
    if kind == 'daemon':
        from . import daemon
        if daemon.status(driver.get('root')) is not None:
            return 'running' if record.get('state') == 'running' else 'queued'
        # no daemon: a queued record waits for the next one, a claimed run is gone
        return 'queued' if record.get('state') == 'queued' else False
    if kind == 'slurm':
        if not driver.get('slurm_id'):
            return False
        from . import scheduler
        try:
            live = scheduler.state(driver['slurm_id'])
        except OSError:
            return None
        return live if live else False
    if kind == 'process':
        if driver.get('host') not in (None, socket.gethostname()):
            return None
        pid = driver.get('pid')
        if pid is None:
            return False
        alive = _process_alive(pid, record.get('id', ''))
        if alive is None:
            return None
        return 'running' if alive else False
    return None


# --- handles ----------------------------------------------------------------------------
class WorkflowHandle:
    """Stable job identity with status, cancellation and durable result retrieval."""

    def __init__(self, folder, future=None, token=None, update=None):
        self.folder = Path(folder).expanduser().resolve()
        self.id = self.folder.name
        self._future = future
        self._token = token
        # writes terminal state for a workflow whose worker never runs (queued cancellation)
        self._update = update
        self._callbacks = []    # done callbacks owed after a queued cancellation

    def _record(self):
        record = _read(self.folder)
        if record is None:
            raise FileNotFoundError(f'No workflow record in {self.folder}')
        if (self.folder / 'result' / 'result.json').is_file():
            record['state'] = 'finished'
        elif record.get('state') in ('queued', 'running') and self._future is None:
            live = driver_state(record)
            if live is False:
                record['state'] = 'interrupted'
            elif live == 'running':
                record['state'] = 'running'
            elif live is None and not record.get('driver'):
                record['state'] = 'interrupted'      # an in-process record of a dead process
        return record

    @property
    def status(self):
        return self._record()['state']

    @property
    def detached(self):
        """The driver kind for a detached workflow, else ``None``."""
        return ((_read(self.folder) or {}).get('driver') or {}).get('kind')

    @property
    def label(self):
        """The display name the submitter gave the workflow, else its id."""
        return self._record().get('label') or self.id

    @property
    def meta(self):
        """The submitter's display annotations (a JSON object), or ``None``."""
        return deepcopy(self._record().get('meta'))

    def _native_runs(self, record):
        """State and progress of each native run this workflow has scheduled so far."""
        from . import jobrecord
        from tumult.launch import monitor
        items = []
        for folder in record.get('runs') or ():
            folder = Path(folder)
            log_path = folder.with_suffix('.log')
            native = jobrecord.read(folder)
            items.append({'output': str(folder),
                          'log': str(log_path) if log_path.is_file() else None,
                          'state': jobrecord.state_of(folder, native),
                          'progress': monitor(folder, log_path if log_path.is_file() else None)})
        return items

    @property
    def progress(self):
        record = self._record()
        diagnostics = deepcopy(record.get('diagnostics') or {})
        phase = record.get('phase')
        if record.get('state') == 'running':
            try:
                native = self._native_runs(record)
            except (OSError, KeyError, ValueError):
                native = []
            if native:
                diagnostics.update(
                    native_runs=len(native),
                    completed=sum(item['state'] == 'finished' for item in native),
                    done=sum(item['progress'][0] for item in native),
                    total=sum(item['progress'][1] for item in native),
                    estimated=True, runs=native)
        return {'phase': phase, 'elapsed': max(
            0.0, (record.get('ended') or time.time()) - record['submitted']),
                'diagnostics': diagnostics}

    @property
    def runs(self):
        return tuple(self._record().get('runs') or ())

    @property
    def logs(self):
        logs = [str(Path(run).with_suffix('.log')) for run in self.runs
                if Path(run).with_suffix('.log').is_file()]
        driver = (self.folder / DRIVER_LOG)
        if driver.is_file():
            logs.append(str(driver))
        return tuple(dict.fromkeys(logs))

    def cancel(self):
        """Cancel this workflow. Returns whether anyone was told.

        In-process: a workflow still queued behind other pool work has its future
        cancelled, so the worker never starts; its record still settles as
        ``cancelled``. A running one is told through its token and settles when the
        worker notices. Detached: a process driver verified to be alive on this host gets
        SIGTERM on its process group and settles its own record; a Slurm driver is
        ``scancel``led, and a job the scheduler never started is settled here, since no
        driver will. An unverified pid is never signalled."""
        if self._future is not None:
            if self._token is None or self._future.done():
                return False
            self._token.cancel()
            if self._future.cancel():
                self._settle_unstarted('cancelled')
            return True
        record = _read(self.folder) or {}
        driver = record.get('driver') or {}
        if record.get('state') not in ('queued', 'running') or not driver:
            return False
        if driver.get('kind') == 'daemon':
            live = driver_state(record)
            if live is False:
                return False
            if live == 'queued' and not (driver.get('pid')):
                from . import daemon
                if daemon.status(driver.get('root')) is None:
                    # nobody will claim it; settle the record here
                    update_record(self.folder, state='cancelled', phase='cancelled',
                                  ended=time.time())
                    return True
            update_record(self.folder, cancel_requested=True)
            return True
        if driver.get('kind') == 'process':
            if driver_state(record) != 'running':
                return False
            os.killpg(int(driver['pid']), signal.SIGTERM)
            self._settle_if_unclaimed()
            return True
        if driver.get('kind') == 'slurm':
            from . import scheduler
            if not driver.get('slurm_id') or not scheduler.scancel(driver['slurm_id']):
                return False
            if _read(self.folder).get('state') == 'queued':
                update_record(self.folder, state='cancelled', phase='cancelled',
                              ended=time.time())
            return True
        return False

    def _settle_if_unclaimed(self):
        """A process driver signalled while its interpreter was still starting dies
        before it can settle anything (the handler is installed after the imports).
        Wait briefly for it to go; if it is gone and the record is still ``queued``,
        nobody else will write, so settle it here. A driver that has claimed the record
        settles it itself."""
        deadline = time.monotonic() + 2 * POLL_S
        while time.monotonic() < deadline:
            time.sleep(0.1)
            record = _read(self.folder) or {}
            if record.get('state') != 'queued':
                return
            if driver_state(record) is False:
                def settle(current):
                    if current.get('state') == 'queued':
                        current.update(state='cancelled', phase='cancelled', ended=time.time())
                mutate_record(self.folder, settle)
                return

    def _settle_unstarted(self, state, error=None):
        """Persist the terminal state of a workflow whose worker never ran, then tell
        its done callbacks; they observe the settled record, not the queued one."""
        values = dict(state=state, phase=state, ended=time.time())
        if error is not None:
            values['error'] = error
        with _lock:
            if self._update is not None:
                self._update(**values)
            _live.pop(self.id, None)
            owed, self._callbacks = self._callbacks, None
        for callback in owed:
            callback(self)

    @property
    def done(self):
        """Whether the workflow has settled, or never will: an in-process future that is
        done, a detached record in a terminal state or without a live driver, or a record
        of a process that no longer exists."""
        if self._future is not None:
            return self._future.done()
        try:
            return self.status in TERMINAL
        except FileNotFoundError:
            return True

    def add_done_callback(self, callback):
        """Call ``callback(handle)`` once the workflow settles, from the worker thread,
        the cancelling one for a queued cancellation, or the watcher for a detached
        workflow; immediately if it already has. Never blocks, and the record is
        terminal by the time the callback runs."""
        if self._future is None:
            if self.done:
                callback(self)
            else:
                _watcher.follow(self, callback)
            return
        with _lock:
            if self._callbacks is not None:
                self._callbacks.append(callback)
                # a worker-run future settles its own record before completing; a
                # cancelled one is reported by _settle_unstarted instead
                self._future.add_done_callback(
                    lambda future: None if future.cancelled() else callback(self))
                return
        callback(self)      # settled by a queued cancellation before registration

    def result(self, timeout=None):
        """The saved result. Raises :class:`CancelledError` for a cancelled workflow and
        ``RuntimeError`` for a failed or interrupted one; waits on a live future, or polls
        a detached record, only if the workflow has not settled, so a done callback may
        call this freely."""
        if self._future is not None:
            try:
                self._future.result(timeout=timeout)
            except TimeoutError:
                raise
            except BaseException:
                pass    # cancelled or failed: the record below names the outcome
        else:
            deadline = None if timeout is None else time.monotonic() + timeout
            while not self.done:
                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(f'Workflow {self.id} is still {self.status}')
                time.sleep(min(POLL_S, 0.2 if timeout is not None and timeout < POLL_S
                               else POLL_S))
        state = self.status
        if state == 'cancelled':
            raise CancelledError(f'Workflow {self.id} was cancelled')
        if state != 'finished':
            raise RuntimeError(f'Workflow {self.id} is {state}: '
                               f'{self._record().get("error") or "no result"}')
        return SynthesisResult.load(self.folder / 'result')

    def to_dict(self):
        return self._record()


class _Watcher:
    """One daemon thread per process that polls followed detached handles and fires
    their callbacks once they settle. Started on the first follow."""

    def __init__(self):
        self._lock = threading.Lock()
        self._followed = []
        self._thread = None

    def follow(self, handle, callback):
        with self._lock:
            self._followed.append((handle, callback))
            if self._thread is None:
                self._thread = threading.Thread(target=self._run, name='tumult-workflow-watch',
                                                daemon=True)
                self._thread.start()

    def _run(self):
        while True:
            time.sleep(POLL_S)
            with self._lock:
                pending = list(self._followed)
            settled = []
            for handle, callback in pending:
                try:
                    if handle.done:
                        settled.append((handle, callback))
                except Exception:
                    log.exception('watching workflow %s', handle.id)
                    settled.append((handle, callback))
            with self._lock:
                for item in settled:
                    self._followed.remove(item)
            for handle, callback in settled:
                try:
                    callback(handle)
                except Exception:
                    log.exception('done callback for workflow %s', handle.id)


_watcher = _Watcher()


# --- execution --------------------------------------------------------------------------
def execute_recorded(folder, resolved, token, base=None):
    """Run ``resolved`` for the record in ``folder``, settling the record as it goes.
    Shared by the in-process worker and the detached driver. Returns the result path;
    raises what the calculation raised after recording it."""
    from .workflows import execute_request
    folder = Path(folder)
    identity = folder.name
    runs = []

    def update(**values):
        return update_record(folder, base, **values)

    def note(output):
        # Called by the batch executor before each native run starts, cache hits included,
        # so the record names its runs while they are in flight.
        if output not in runs:
            runs.append(output)
            update(runs=list(runs))

    update(state='running', phase='preparation', started=time.time())
    try:
        token.check()
        update(phase='synthesis')
        with observe_runs(note):
            result = execute_request(resolved, cancel_token=token)
        update(phase='saving')
        result_path = result.save(folder / 'result')
        update(state='finished', phase='finished', ended=time.time(),
               result=str(result_path), runs=list(dict.fromkeys(runs + list(result.runs))))
        return str(result_path)
    except CancelledError:
        update(state='cancelled', phase='cancelled', ended=time.time())
        raise
    except BaseException as exc:
        # the traceback belongs to the process log whichever document, if any, is
        # still waiting; the record keeps the concise cause for its Processes row
        log.exception('workflow %s failed', identity)
        update(state='failed', phase='failed', ended=time.time(),
               error=f'{type(exc).__name__}: {exc}')
        raise


def _start_driver(folder, record, detach, slurm_script):
    """Start the detached driver for ``folder`` and record how to find it again."""
    if os.name != 'posix':
        raise ValueError('Detached workflows need a POSIX host')
    env = {**os.environ}
    package_root = str(Path(__file__).resolve().parents[2])
    env['PYTHONPATH'] = os.pathsep.join(
        part for part in (package_root, env.get('PYTHONPATH')) if part)
    driver = {'kind': detach, 'pid': None, 'host': None, 'slurm_id': None,
              'script': None, 'cwd': os.getcwd(), 'python': sys.executable}
    if detach == 'daemon':
        from . import daemon
        driver['root'] = str(folder.parent.parent)
        update_record(folder, record, driver=driver)     # visible before the daemon looks
        daemon.ensure(driver['root'])
        return _read(folder)
    if detach == 'process':
        log_path = folder / DRIVER_LOG
        with open(log_path, 'ab') as log_file:
            child = subprocess.Popen(
                [sys.executable, '-m', 'tumult.synthesis.driver', str(folder)],
                cwd=driver['cwd'], env=env, stdin=subprocess.DEVNULL,
                stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True)
        driver.update(pid=child.pid, host=socket.gethostname())
    else:
        if not slurm_script:
            raise ValueError("detach='slurm' needs the sbatch script that runs the driver")
        from . import scheduler
        script = Path(slurm_script).expanduser()
        if not script.is_file():
            raise FileNotFoundError(f'Slurm script not found: {script}')
        driver.update(script=str(script.resolve()),
                      slurm_id=scheduler.sbatch(script, [str(folder)], env=env))
    return update_record(folder, record, driver=driver)


def submit(request, *, root=None, label=None, meta=None, detach=None, slurm_script=None):
    """Copy, resolve and record ``request``; return immediately with a durable handle.

    ``label`` and ``meta`` are the submitter's display name and annotations (JSON), kept
    in the record so a later process can list and present the workflow as it was
    submitted. ``detach`` picks the driver: ``None`` runs in this process, ``'daemon'``
    hands it to the root's worker daemon (started if needed), ``'process'`` starts a
    driver that outlives it, ``'slurm'`` submits that driver through
    ``sbatch slurm_script <folder>``. Absent an argument, ``request.prefs`` may name the
    driver as ``workflow_runner`` (``local`` means the daemon, ``in-process`` none) and
    the script as ``slurm_script``."""
    from .workflows import resolve_request
    resolved = request if isinstance(request, ResolvedSynthRequest) else resolve_request(request)
    # JSON round-trip proves that no mutable masks/settings remain shared with the caller.
    owned = SynthRequest.from_dict(resolved.request.to_dict())
    resolved = resolve_request(owned)
    prefs = owned.prefs or {}
    if detach is None:
        choice = str(prefs.get('workflow_runner') or '')
        detach = RUNNER_PREFS.get(choice, choice)
    if detach is not None and detach not in DETACHED:
        raise ValueError(f'detach must be one of {DETACHED} or None, not {detach!r}')
    if detach == 'slurm' and slurm_script is None and prefs.get('slurm_script'):
        # a relative script path in the preferences is relative to the dispatch directory
        slurm_script = Path(prefs.get('dispatch_cwd') or '.').expanduser() / prefs['slurm_script']
    base = Path(root or resolved.output_root) / '_workflows'
    identity = f'{time.strftime("%Y%m%dT%H%M%S")}-{uuid.uuid4().hex[:10]}'
    folder = base / identity
    record = {'version': 1, 'id': identity, 'state': 'queued', 'phase': 'queued',
              'submitted': time.time(), 'request': resolved.to_dict(), 'runs': [],
              'result': None, 'label': label, 'meta': deepcopy(meta), 'driver': None}
    _atomic(folder / RECORD, record)

    if detach is not None:
        try:
            _start_driver(folder, record, detach, slurm_script)
        except BaseException as exc:
            update_record(folder, record, state='failed', phase='failed', ended=time.time(),
                          error=f'{type(exc).__name__}: {exc}')
            raise
        return WorkflowHandle(folder)

    token = CancelToken()

    def update(**values):
        return update_record(folder, record, **values)

    handle = WorkflowHandle(folder, token=token, update=update)

    def work():
        try:
            return execute_recorded(folder, resolved, token, record)
        finally:
            with _lock:
                _live.pop(identity, None)

    with _lock:
        _live[identity] = handle
        try:
            handle._future = WORKFLOWS.submit(work)
        except BaseException as exc:
            handle._settle_unstarted('failed', error=f'{type(exc).__name__}: {exc}')
            raise
    return handle


def cancel_live_workflows():
    """Cancel every in-process workflow this process still holds live, for server
    shutdown; returns how many were told. Queued ones settle as ``cancelled`` at once,
    running ones when their worker notices the token. Detached workflows are not held
    here and keep running."""
    with _lock:
        handles = list(_live.values())
    return sum(handle.cancel() for handle in handles)


def open_workflow(path, *, root=None):
    """Open a workflow by folder/path or stable id after a UI reconnect."""
    path = Path(path)
    folder = ((Path(root) / '_workflows' / path) if root is not None and len(path.parts) == 1
              else path).expanduser().resolve()
    with _lock:
        return _live.get(folder.name) or WorkflowHandle(folder)


def discover_workflows(root, limit=None, revive=True):
    """Newest workflow handles under ``<root>/_workflows``, whatever their state: live
    ones (this process holds their handle, or their detached driver still runs),
    finished results, and records left queued or running by a driver that no longer
    exists, which report ``interrupted``. Only that folder is scanned, never the run
    cache beside it. With ``revive``, a queued daemon workflow whose daemon has gone
    gets a new daemon started."""
    base = Path(root).expanduser().resolve() / '_workflows'
    folders = [path.parent for path in base.glob(f'*/{RECORD}')]
    folders.sort(key=lambda path: (_read(path) or {}).get('submitted', 0), reverse=True)
    handles = [open_workflow(path) for path in folders]
    if revive and any(_waits_for_daemon(handle) for handle in handles):
        from . import daemon
        daemon.ensure(base.parent)      # queued daemon work with no daemon: start one
    return handles[:limit] if limit else handles


def _waits_for_daemon(handle):
    record = _read(handle.folder) or {}
    driver = record.get('driver') or {}
    return (driver.get('kind') == 'daemon' and record.get('state') == 'queued'
            and not driver.get('pid'))


def run(request, *, cancel_token=None):
    """Blocking counterpart to :func:`submit`."""
    from .workflows import execute_request
    return execute_request(request, cancel_token=cancel_token)
