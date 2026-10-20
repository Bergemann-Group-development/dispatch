"""Durable asynchronous submission of the same calculations used synchronously.

A workflow is recorded in ``<root>/_workflows/<id>/workflow.json``: one public operation
of :mod:`tumult.synthesis.tasks` on a request (:func:`submit_operation`, the record's
``task`` block). The submitter writes the record ``queued`` and hands it to one of two
drivers:

* ``daemon`` (the default): the workflow daemon of this user and host
  (:mod:`tumult.synthesis.daemon`), started on demand. It runs workflows one at a time
  with the whole CPU budget, outlives the submitter, and reports state, threads and
  progress to every handle as events.
* ``slurm``: ``python -m tumult.synthesis.driver <folder>`` submitted through
  ``sbatch <script> <folder>``. The scheduler owns liveness (``squeue``) and cancellation
  (``scancel``); a handle polls the record.

The driver settles the record, and every update takes the folder's file lock. Scripts
that want no daemon call the operation itself instead: blocking, in-process and not
recorded as a workflow.
"""
from __future__ import annotations

from copy import deepcopy
from contextlib import contextmanager
import json
import logging
import os
from pathlib import Path
import socket
import sys
import threading
import time
import uuid

from .runtime import observe_runs, workflow_scope
from .cancel import CancelledError
from . import tasks

RECORD = 'workflow.json'
LOCK = 'workflow.lock'
DRIVER_LOG = 'driver.log'
DETACHED = ('daemon', 'slurm')
# preference spellings of the driver choice (TumultStudio2 Preferences -> Engine)
RUNNER_PREFS = {'local': 'daemon', 'daemon': 'daemon', 'slurm': 'slurm', '': 'daemon'}
TERMINAL = ('finished', 'failed', 'cancelled', 'interrupted')
POLL_S = 3.0        # how often a Slurm record is re-read, and the fallback wait for results
_lock = threading.RLock()
log = logging.getLogger(__name__)


# --- the record on disk ----------------------------------------------------------------
def _atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'{path.name}.{os.getpid()}.{threading.get_ident()}.partial')
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
    file system) serialises with the submitter."""
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
    """The live state of a record's driver: ``queued``, ``running``, ``False`` (gone) or
    ``None`` (not observable from here: a daemon on another host, or no scheduler
    client). A daemon record waiting in the queue stays ``queued`` whether or not a
    daemon runs now; the next daemon to register its root runs it."""
    driver = record.get('driver') or {}
    kind = driver.get('kind')
    if kind == 'daemon':
        if record.get('state') == 'queued':
            return 'queued'
        from . import daemon
        live = daemon.runs_record(driver)
        return None if live is None else ('running' if live else False)
    if kind == 'slurm':
        if not driver.get('slurm_id'):
            return False
        from . import scheduler
        try:
            live = scheduler.state(driver['slurm_id'])
        except OSError:
            return None
        return live if live else False
    return False        # no driver this code knows: the process that ran it is gone


# --- handles ----------------------------------------------------------------------------
class WorkflowHandle:
    """Stable job identity with status, cancellation and durable result retrieval."""

    def __init__(self, folder):
        self.folder = Path(folder).expanduser().resolve()
        self.id = self.folder.name

    def _record(self):
        record = _read(self.folder)
        if record is None:
            raise FileNotFoundError(f'No workflow record in {self.folder}')
        if (self.folder / 'result' / 'result.json').is_file():
            record['state'] = 'finished'
        elif record.get('state') in ('queued', 'running'):
            live = driver_state(record)
            if live is False:
                record['state'] = 'interrupted'
            elif live == 'running':
                record['state'] = 'running'
        return record

    @property
    def status(self):
        return self._record()['state']

    @property
    def detached(self):
        """The driver kind: ``daemon`` or ``slurm``."""
        return ((_read(self.folder) or {}).get('driver') or {}).get('kind')

    @property
    def kind(self):
        """The :class:`tasks.TaskKind` of the calculation this record carries."""
        return tasks.kind_of(self._record())

    def progress_stamp(self):
        """A value that changes whenever a new progress snapshot is written (its mtime),
        or ``None`` before the first."""
        try:
            return (self.folder / tasks.PROGRESS_FILE).stat().st_mtime_ns
        except FileNotFoundError:
            return None

    def progress_result(self):
        """The latest intermediate result the calculation reported, decoded, or ``None``
        for a kind without progress or before the first report."""
        kind = self.kind
        return tasks.read_progress(self.folder, kind) if kind.progress else None

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
        from tumult.launch import monitor, run_file
        items = []
        for folder in record.get('runs') or ():
            folder = Path(folder)
            log_path = run_file(folder, '.log')
            native = jobrecord.read(folder)
            items.append({'output': str(folder),
                          'log': str(log_path) if log_path.is_file() else None,
                          'state': jobrecord.state_of(folder, native),
                          'progress': monitor(folder, log_path if log_path.is_file() else None)})
        return items

    def native_fraction(self, since=None):
        """How far the dispatch runs of the step under way are, as the mean over that
        step's runs of each one's fraction: a run still to start counts 0, a running one
        its ``done / total`` from :func:`tumult.launch.monitor` (an LTE synthesis counts
        its RT patch tasks, a departure solve its iterations), a finished one 1. A step's
        runs are the ones the record names that have not finished before ``since`` (by
        default the last progress snapshot, else the record's start): the batch executor
        names every run of a batch before launching any, so the mean only grows within a
        step. A run finished earlier is a cache hit and costs nothing; a failed or
        cancelled one is left out. ``None`` while the step has no runs yet. Cheap enough
        to poll every few seconds: the record is read as written (no driver check), only
        the launcher record of each run is read, and the queue and dispatch logs are
        small text files."""
        from . import jobrecord
        from tumult.launch import monitor, run_file
        record = _read(self.folder) or {}
        if since is None:
            stamp = self.progress_stamp()
            since = stamp / 1e9 if stamp else float(record.get('started') or 0.0)
        fractions = []
        for folder in record.get('runs') or ():
            folder = Path(folder)
            launcher = jobrecord.read_launcher(folder)
            if launcher is None:
                fractions.append(0.0)                 # named, not launched yet
                continue
            state = launcher.get('state')
            if state == 'finished':
                if float(launcher.get('ended') or 0.0) >= since:
                    fractions.append(1.0)
                continue
            if state != 'running':
                continue
            log = run_file(folder, '.log')
            done, total = monitor(folder, log if log.is_file() else None)
            fractions.append(min(1.0, float(done) / float(total)) if total else 0.0)
        return sum(fractions) / len(fractions) if fractions else None

    @property
    def progress(self):
        record = self._record()
        diagnostics = deepcopy(record.get('diagnostics') or {})
        phase = record.get('phase')
        native = []
        if record.get('state') == 'running':
            try:
                native = self._native_runs(record)
            except (OSError, KeyError, ValueError):
                pass
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
        from tumult.launch import run_file
        logs = [str(run_file(run, '.log')) for run in self.runs
                if run_file(run, '.log').is_file()]
        driver = (self.folder / DRIVER_LOG)
        if driver.is_file():
            logs.append(str(driver))
        return tuple(dict.fromkeys(logs))

    def cancel(self):
        """Cancel this workflow; returns whether anyone was told. A daemon workflow is
        cancelled by the daemon: a queued one settles at once, a running one once its
        worker notices the token. A Slurm job is ``scancel``led, and one the scheduler
        never started is settled here, since no driver will. A record whose driver is
        gone is left alone."""
        try:
            record = self._record()
        except FileNotFoundError:
            return False
        driver = record.get('driver') or {}
        if record['state'] not in ('queued', 'running'):
            return False
        if driver.get('kind') == 'daemon':
            from . import daemon
            reply = daemon.client().request('cancel', job=self.id, folder=str(self.folder))
            return reply.get('state') in ('cancelled', 'cancelling')
        if driver.get('kind') == 'slurm':
            from . import scheduler
            if not driver.get('slurm_id') or not scheduler.scancel(driver['slurm_id']):
                return False
            if (_read(self.folder) or {}).get('state') == 'queued':
                update_record(self.folder, state='cancelled', phase='cancelled',
                              ended=time.time())
            return True
        return False

    @property
    def done(self):
        """Whether the workflow has settled, or never will (its driver is gone)."""
        try:
            return self.status in TERMINAL
        except FileNotFoundError:
            return True

    def add_done_callback(self, callback):
        """Call ``callback(handle)`` once the workflow settles: at once if it has, else
        from the daemon client's event thread (a daemon workflow) or the polling watcher
        (a Slurm job). Never blocks, and the record is terminal when the callback runs."""
        if self.done:
            callback(self)
            return
        if self.detached == 'slurm':
            _watcher.follow(self, callback)
            return
        from . import daemon
        client = daemon.client()
        fired, unfollow = [False], [None]

        def check(event=None):
            if event is not None and event.get('event') not in ('state', 'reconnect'):
                return
            if not self.done:
                return
            with _lock:
                if fired[0]:
                    return
                fired[0] = True
            if unfollow[0] is not None:
                unfollow[0]()
            callback(self)

        unfollow[0] = client.follow(self.id, check)
        if fired[0]:
            unfollow[0]()
        if (_read(self.folder) or {}).get('state') == 'queued':
            # a record queued while no daemon ran: make sure this daemon knows it
            client.request('submit', folder=str(self.folder))
        check()

    def on_progress(self, callback):
        """Call ``callback(stamp)`` whenever the daemon reports a new progress snapshot
        (and once after a reconnect, which may have missed one). Returns a function that
        stops it. Slurm records report nothing here; poll :meth:`progress_stamp`."""
        if self.detached != 'daemon':
            return lambda: None
        from . import daemon

        def relay(event):
            if event.get('event') == 'progress':
                callback(event.get('stamp'))
            elif event.get('event') == 'reconnect':
                callback(self.progress_stamp())

        return daemon.client().follow(self.id, relay)

    def result(self, timeout=None):
        """The saved result. Raises :class:`CancelledError` for a cancelled workflow and
        ``RuntimeError`` for a failed or interrupted one; waits for the workflow to settle
        first (``TimeoutError`` after ``timeout`` seconds), so a done callback may call
        this freely."""
        if not self.done:
            settled = threading.Event()
            self.add_done_callback(lambda handle: settled.set())
            deadline = None if timeout is None else time.monotonic() + timeout
            while not self.done:
                remaining = POLL_S if deadline is None else deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f'Workflow {self.id} is still {self.status}')
                settled.wait(min(POLL_S, remaining))
        state = self.status
        if state == 'cancelled':
            raise CancelledError(f'Workflow {self.id} was cancelled')
        if state != 'finished':
            raise RuntimeError(f'Workflow {self.id} is {state}: '
                               f'{self._record().get("error") or "no result"}')
        return self.kind.load(self.folder / tasks.RESULT_DIR)

    def to_dict(self):
        return self._record()


class _Watcher:
    """One thread per process that polls followed Slurm handles and fires their
    callbacks once they settle. Started on the first follow."""

    def __init__(self):
        self._lock = threading.Lock()
        self._followed = []
        self._thread = None

    def follow(self, handle, callback):
        with self._lock:
            self._followed.append((handle, callback))
            if self._thread is None:
                self._thread = threading.Thread(target=self._run, name='tumult-slurm-watch',
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
def decode_task(record):
    """``(kind, args)`` of a record: the task kind and its arguments decoded into the
    live objects ``kind.run`` takes. Raises what the decoding raised."""
    kind = tasks.kind_of(record)
    return kind, kind.decode(tasks.args_of(record))


PROGRESS_MIN_S = 1.0    # how often at most a running calculation rewrites progress.json


def execute_recorded(folder, kind, args, token, base=None, listener=None):
    """Run ``kind`` on ``args`` for the record in ``folder``, settling the record as it
    goes. Shared by the daemon and the Slurm driver. Intermediate results the
    calculation reports are written as the record's progress snapshot, at most every
    ``PROGRESS_MIN_S``. The leaf jobs it submits are labelled with the workflow id (see
    :func:`runtime.workflow_scope`). ``listener(record)`` sees every record update as
    written and ``listener(None, progress_stamp=...)`` every progress snapshot. Returns
    the result path; raises what the calculation raised after recording it."""
    folder = Path(folder)
    identity = folder.name
    runs = []
    last_progress = [0.0]

    def update(**values):
        record = update_record(folder, base, **values)
        if listener is not None:
            listener(record)
        return record

    def note(output):
        # Called by the batch executor before each native run starts, cache hits included,
        # so the record names its runs while they are in flight.
        if output not in runs:
            runs.append(output)
            update(runs=list(runs))

    def progress(result):
        now = time.monotonic()
        if now - last_progress[0] < PROGRESS_MIN_S:
            return
        last_progress[0] = now
        try:
            tasks.write_progress(folder, kind, result)
        except OSError:
            log.exception('workflow %s: progress snapshot not written', identity)
            return
        if listener is not None:
            try:
                stamp = (folder / tasks.PROGRESS_FILE).stat().st_mtime_ns
            except OSError:
                return
            listener(None, progress_stamp=stamp)

    update(state='running', phase='preparation', started=time.time())
    try:
        token.check()
        update(phase='synthesis')
        with observe_runs(note), workflow_scope(identity):
            result = kind.run(args, token, progress)
        update(phase='saving')
        result_path = kind.save(result, folder / tasks.RESULT_DIR)
        update(state='finished', phase='finished', ended=time.time(),
               result=str(result_path),
               runs=list(dict.fromkeys(runs + list(getattr(result, 'runs', ())))))
        return str(result_path)
    except CancelledError:
        update(state='cancelled', phase='cancelled', ended=time.time())
        raise
    except BaseException as exc:
        log.exception('workflow %s failed', identity)
        update(state='failed', phase='failed', ended=time.time(),
               error=f'{type(exc).__name__}: {exc}')
        raise


def _submit_slurm(folder, record, slurm_script):
    """Submit the Slurm driver for ``folder`` and record its job id."""
    if not slurm_script:
        raise ValueError("detach='slurm' needs the sbatch script that runs the driver")
    from . import scheduler
    script = Path(slurm_script).expanduser()
    if not script.is_file():
        raise FileNotFoundError(f'Slurm script not found: {script}')
    env = {**os.environ}
    package_root = str(Path(__file__).resolve().parents[2])
    env['PYTHONPATH'] = os.pathsep.join(
        part for part in (package_root, env.get('PYTHONPATH')) if part)
    driver = {'kind': 'slurm', 'pid': None, 'host': None, 'script': str(script.resolve()),
              'cwd': os.getcwd(), 'python': sys.executable,
              'slurm_id': scheduler.sbatch(script, [str(folder)], env=env)}
    return update_record(folder, record, driver=driver)


def _driver_choice(prefs, detach, slurm_script):
    """The driver a submission uses: the explicit ``detach``, else the preferences'
    ``workflow_runner``; the Slurm script likewise falls back to ``slurm_script``."""
    prefs = prefs or {}
    if detach is None:
        choice = str(prefs.get('workflow_runner') or '')
        if choice not in RUNNER_PREFS:
            raise ValueError(f'workflow_runner must be one of {sorted(RUNNER_PREFS)}, '
                             f'not {choice!r}')
        detach = RUNNER_PREFS[choice]
    if detach not in DETACHED:
        raise ValueError(f'detach must be one of {DETACHED} or None, not {detach!r}')
    if detach == 'slurm' and slurm_script is None and prefs.get('slurm_script'):
        # a relative script path in the preferences is relative to the dispatch directory
        slurm_script = Path(prefs.get('dispatch_cwd') or '.').expanduser() / prefs['slurm_script']
    return detach, slurm_script


def _submit_record(base, record, detach, slurm_script):
    """Write ``record`` under ``base`` and hand it to its driver; returns the handle. A
    submission that fails leaves the record ``failed`` with the cause."""
    folder = Path(base).expanduser().resolve() / record['id']
    if detach == 'daemon':
        record['driver'] = {'kind': 'daemon', 'host': socket.gethostname(), 'pid': None}
    _atomic(folder / RECORD, record)
    try:
        if detach == 'slurm':
            _submit_slurm(folder, record, slurm_script)
        else:
            from . import daemon
            daemon.client().request('submit', folder=str(folder))
    except BaseException as exc:
        update_record(folder, record, state='failed', phase='failed', ended=time.time(),
                      error=f'{type(exc).__name__}: {exc}')
        raise
    return WorkflowHandle(folder)


def _new_record(label, meta):
    identity = f'{time.strftime("%Y%m%dT%H%M%S")}-{uuid.uuid4().hex[:10]}'
    return {'version': 1, 'id': identity, 'state': 'queued', 'phase': 'queued',
            'submitted': time.time(), 'task': None, 'runs': [],
            'result': None, 'label': label, 'meta': deepcopy(meta), 'driver': None}


def submit_operation(operation, request, *, label=None, meta=None, root=None, detach=None,
                     slurm_script=None, **kwargs):
    """Record a public operation (``"synthesis.synthesize"``, ``"synthesis.fit"``,
    ``"departures.solve"``, ``"abundance.fit"``, ``"abundance.nlte_corrections"`` or
    ``"sparams.ionbalance"``) on ``request`` with its keyword arguments, and return a
    durable handle at once.

    The driver runs the same function a script calls, with cancellation and progress;
    ``handle.result()`` returns what that call returns. ``label`` and ``meta`` are the
    submitter's display name (the operation by default) and annotations (JSON), kept in
    the record so a later process can list and present the workflow as it was submitted.
    The record lives under ``root``, by default the request's cache
    (``<config.root>/<config.cache>``). ``detach`` picks the driver: ``'daemon'`` (the
    default) hands it to this host's workflow daemon, ``'slurm'`` submits the driver
    through ``sbatch slurm_script <folder>``; absent an argument, the configuration's
    preferences may name the driver as ``workflow_runner`` (``local`` means the daemon)
    and the script as ``slurm_script``. The arguments are decoded here first, so an
    unreadable argument fails the caller rather than the driver."""
    from .workflows import output_root
    task = tasks.KINDS.get(operation)
    if task is None:
        raise ValueError(f'unknown operation {operation!r}; one of {sorted(tasks.KINDS)}')
    prefs = request.config.to_prefs()
    args = json.loads(json.dumps(tasks.Operation.encode_args(request, kwargs),
                                 default=tasks._json_default))   # owned, JSON-clean
    task.decode(args)
    detach, slurm_script = _driver_choice(prefs, detach, slurm_script)
    record = _new_record(label or operation, meta)
    record['task'] = {'kind': task.name, 'args': args}
    base = Path(root or output_root(prefs)) / '_workflows'
    return _submit_record(base, record, detach, slurm_script)


def open_workflow(path, *, root=None):
    """Open a workflow by folder/path, or by stable id under ``root``."""
    path = Path(path)
    folder = ((Path(root) / '_workflows' / path) if root is not None and len(path.parts) == 1
              else path)
    return WorkflowHandle(folder)


def discover_workflows(root, limit=None):
    """Newest workflow handles under ``<root>/_workflows``, whatever their state: running
    and queued ones, finished results, and records whose driver is gone, which report
    ``interrupted``. Only that folder is scanned, never the run cache beside it."""
    base = Path(root).expanduser().resolve() / '_workflows'
    folders = [path.parent for path in base.glob(f'*/{RECORD}')]
    folders.sort(key=lambda path: (_read(path) or {}).get('submitted', 0), reverse=True)
    handles = [WorkflowHandle(path) for path in folders]
    return handles[:limit] if limit else handles
