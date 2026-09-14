"""One workflow daemon per user and host, the only executor of durable workflows.

    python -m tumult.synthesis.daemon [max_threads]

Bokeh servers and scripts connect to it over a Unix socket in ``/tmp/tumult-<uid>/``
(mode 0700; ``$TUMULT_DAEMON_DIR`` names another directory). A per-user directory under
``/tmp`` keeps the socket off NFS home directories, which do not carry Unix sockets, and
within the 104-byte socket path limit that macOS's ``$TMPDIR`` exceeds. The daemon owns
the CPU budget (:class:`runtime.ExecutionService`), runs one workflow at a time in FIFO
order with the full budget, and reports per-workflow threads, the budget in use and the
queue to every connection, so a job reaches its submitter as soon as it settles.

Protocol: newline-delimited JSON. A request carries an ``id`` that its reply echoes; a
failed request is answered with ``error``. Events carry ``event`` and no ``id``.

========================================  ===================================================
request                                   reply
========================================  ===================================================
``hello {role, pid, host, roots,``        ``{daemon_pid, host, budget}``; first message.
``max_threads}``                          ``roots`` are registered for rescans; a server's
                                          ``max_threads`` configures the budget.
``submit {folder}``                       ``{accepted, state}``; the client wrote the record
                                          ``queued`` first, so the job is durable already.
``cancel {job, folder}``                  ``{state}``; a workflow id or an external run id.
``list {roots, limit}``                   ``{jobs, budget, used}``: running, queued, external
                                          runs and the ``limit`` newest settled records of
                                          every registered root; ``roots`` are registered.
``subscribe {jobs}``                      ``{}``; ``jobs`` is a list of ids or ``"all"``.
``configure {max_threads}``               ``{budget}``; lowering the budget drains.
========================================  ===================================================

Events: ``state`` (a job summary: ``job``, ``state``, ``phase``, ``error``, ``label``,
``kind``, ``folder`` and the timestamps), ``threads {job, threads}``, ``progress {job,
stamp}``, ``usage {used, budget}`` and ``external {jobs}`` when the live external runs
change. A progress event carries only the stamp; the snapshot itself stays in
``progress.json``, since it can hold a spectrum.

Records stay under each job's output root and remain the journal. After ``submit`` the
daemon is their only writer. Registering a root re-queues its ``queued`` daemon records
and marks ``running`` ones this daemon does not run ``interrupted``: the daemon that ran
them has exited. The registered roots are kept in ``roots.json`` for the next daemon.

Lifetime: the daemon exits once it has no running or queued workflow and no connection
with the ``server`` role. A ``client`` connection (a script, a notebook) does not keep
it alive, but the work it submitted does, so a script may submit and exit. A daemon that
never had a connection waits ``GRACE_S`` first, so its spawner's ``hello`` cannot lose
the race. To exit it takes ``start.lock``, closes and unlinks the socket, re-checks for
work or a server that arrived in between (and listens again if so), and releases
``daemon.lock`` before ``start.lock``; a client that cannot connect after that starts a
new daemon. SIGTERM, SIGINT and SIGHUP cancel the running workflow and exit; queued
records wait for the next daemon.

While a server is connected, the daemon scans the registered roots every ``SCAN_S`` for
dispatch runs recorded as live (``job.nml``/``job.json``, see :mod:`jobrecord`) that no
workflow record names: runs started from a shell, a notebook or
:func:`submission.run`. They are listed as ``external`` jobs and cancelled through
:func:`jobrecord.cancel`.
"""
from __future__ import annotations

from contextlib import contextmanager
import itertools
import json
import logging
import os
from pathlib import Path
import queue
import select
import signal
import socket
import subprocess
import sys
import threading
import time

from .cancel import CancelToken, CancelledError
from .runtime import WORKFLOWS, get_runtime, resolve_cpu_budget
from . import jobrecord, submission, tasks

log = logging.getLogger(__name__)

SOCKET = 'daemon.sock'
LOCKFILE = 'daemon.lock'        # held for the daemon's lifetime
START_LOCK = 'start.lock'       # serialises starting a daemon with a daemon exiting
INFO = 'daemon.json'
ROOTS = 'roots.json'
LOGFILE = 'daemon.log'
GRACE_S = 30.0          # how long a daemon nobody connected to yet waits before exiting
CHECK_S = 5.0           # how often the exit rule is checked without another reason to
RECENT_S = 0.5          # how long after a connection or request the daemon stays anyway
SCAN_S = 10.0           # external-run scan period while a server is connected
START_TIMEOUT_S = 10.0  # how long a client waits for a daemon it started
LIST_LIMIT = 20         # settled records per root in a ``list`` reply
TERMINAL = submission.TERMINAL
ROLES = ('server', 'client')


def directory():
    """The per-user directory holding the socket, the locks and the log, created with
    mode 0700. Refuses a directory another user owns or others can enter."""
    configured = os.environ.get('TUMULT_DAEMON_DIR')
    path = Path(configured) if configured else Path('/tmp') / f'tumult-{os.getuid()}'
    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    status = path.stat()
    if status.st_uid != os.getuid() or status.st_mode & 0o077:
        raise PermissionError(f'{path} must belong to this user with mode 0700')
    return path


@contextmanager
def _flocked(path, shared=False, block=True):
    """Hold ``flock`` on ``path``; raises ``BlockingIOError`` without ``block`` when it
    is taken."""
    import fcntl
    with open(path, 'a+') as handle:
        mode = (fcntl.LOCK_SH if shared else fcntl.LOCK_EX) | (0 if block else fcntl.LOCK_NB)
        fcntl.flock(handle, mode)
        try:
            yield handle
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def running():
    """The live daemon's ``daemon.json`` (``pid``, ``host``, ``budget``), else ``None``.
    Liveness is the daemon's lock, never a pid, so a reused pid cannot pass for it."""
    base = directory()
    try:
        with _flocked(base / LOCKFILE, shared=True, block=False):
            return None
    except BlockingIOError:
        pass
    try:
        return json.loads((base / INFO).read_text())
    except (OSError, ValueError):
        return {}           # starting: the lock is held, the info not yet written


def runs_record(driver):
    """Whether the daemon that claimed a record (its ``driver`` block) is still running:
    ``True``/``False`` on this host, ``None`` for a record claimed on another host."""
    if driver.get('host') not in (None, socket.gethostname()):
        return None
    info = running()
    return bool(info) and info.get('pid') == driver.get('pid')


def _summary(folder, record, threads=0):
    """The job description used by ``list`` and ``state`` events."""
    driver = record.get('driver') or {}
    try:
        kind = tasks.kind_of(record).name
    except KeyError:
        kind = (record.get('task') or {}).get('kind')
    return {'job': record.get('id') or Path(folder).name, 'type': 'workflow',
            'folder': str(folder), 'label': record.get('label'), 'kind': kind,
            'state': record.get('state'), 'phase': record.get('phase'),
            'error': record.get('error'), 'submitted': record.get('submitted'),
            'started': record.get('started'), 'ended': record.get('ended'),
            'driver': driver.get('kind'), 'slurm_id': driver.get('slurm_id'),
            'threads': threads, 'meta': record.get('meta')}


def _external_summary(record):
    folder = Path(record['output'])
    threads = None
    if record.get('host') in (None, socket.gethostname()) and not record.get('slurm_id'):
        threads = int(record.get('omp_threads') or 1) * int(record.get('mpi_size') or 1)
    return {'job': jobrecord.run_id(folder), 'type': 'external', 'folder': str(folder),
            'label': folder.name, 'kind': 'dispatch', 'state': record['resolved'],
            'phase': None, 'error': None, 'started': jobrecord.started_time(record),
            'pid': record.get('pid'), 'host': record.get('host'),
            'slurm_id': record.get('slurm_id'), 'input': record.get('input'),
            'threads': threads}


# --- the daemon ---------------------------------------------------------------------------
class _Connection:
    """One accepted connection: a reader thread handling requests and a writer thread
    draining an outbound queue, so a slow peer never blocks the daemon."""
    QUEUE = 10000

    def __init__(self, sock):
        self.sock = sock
        self.role = None
        self.pid = None
        self.all = False
        self.subscribed = set()
        self.accepted = time.monotonic()
        self.outbox = queue.Queue(self.QUEUE)
        self.closed = False

    def send(self, message):
        if self.closed:
            return
        try:
            self.outbox.put_nowait(message)
        except queue.Full:
            log.warning('connection of pid %s is not reading; closing it', self.pid)
            self.close()

    def wants(self, job):
        return self.all or job in self.subscribed

    def write_loop(self):
        while True:
            message = self.outbox.get()
            if message is None:
                return
            try:
                self.sock.sendall((json.dumps(message, default=str) + '\n').encode())
            except OSError:
                self.close()
                return

    def close(self):
        if self.closed:
            return
        self.closed = True
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except OSError:
            pass
        try:
            self.outbox.put_nowait(None)
        except queue.Full:
            pass


class _Workflow:
    def __init__(self, folder, record):
        self.folder = Path(folder)
        self.id = record.get('id') or self.folder.name
        self.submitted = record.get('submitted') or 0.0
        self.record = record
        self.token = None
        self.future = None
        self.threads = 0
        self.emitted = None     # the (state, phase) last sent in a state event


class Daemon:
    def __init__(self, max_threads=0, *, grace_s=GRACE_S, exit_when_idle=True):
        self.budget = resolve_cpu_budget(max_threads)
        self.grace_s = grace_s
        self.exit_when_idle = exit_when_idle
        self._lock = threading.RLock()      # queue, active workflow, roots
        self._conns_lock = threading.Lock()  # the connection list only; never calls out
        self._connections = []
        self._queue = []                    # _Workflow, FIFO by (submitted, id)
        self._active = None
        self._roots = set()
        self._external = {}                 # run id -> summary
        self._external_cache = {}           # record path -> (mtime_ns, record)
        self._stop = threading.Event()
        self._wake = threading.Event()
        self._ever_connected = False
        self._started = time.monotonic()
        self._last_request = 0.0
        self._listener = None
        self._threads = []

    # --- lifetime ------------------------------------------------------------------
    def serve(self):
        """Run until the exit rule holds or a signal arrives; returns the exit status."""
        import fcntl
        self.dir = directory()
        self._lockfile = open(self.dir / LOCKFILE, 'a+')
        deadline = time.monotonic() + 2.0   # a liveness check holds the lock shared briefly
        while True:
            try:
                fcntl.flock(self._lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() > deadline or running():
                    log.info('another daemon holds %s', self.dir / LOCKFILE)
                    self._lockfile.close()
                    return 0
                time.sleep(0.05)
        start_lock = None
        try:
            self._listen()
            if threading.current_thread() is threading.main_thread():
                for name in ('SIGTERM', 'SIGINT', 'SIGHUP'):
                    if hasattr(signal, name):
                        signal.signal(getattr(signal, name), self._on_signal)
            runtime = get_runtime()
            runtime.configure(self.budget)
            runtime.set_listener(self._on_allocation)
            self._write_info()
            for root in self._stored_roots():
                self._register(root)
            self._spawn(self._scan_loop, 'tumult-daemon-scan')
            log.info('serving %s as pid %d with budget %d', self.dir / SOCKET, os.getpid(),
                     self.budget)
            wait = CHECK_S
            while not self._stop.is_set():
                self._wake.wait(wait)
                self._wake.clear()
                wait = self._idle() if self.exit_when_idle else CHECK_S
                if wait == 0:
                    start_lock = self._close_if_idle()
                    if start_lock is not None:
                        log.info('no work and no server, exiting')
                        break
                    wait = CHECK_S
            if self._stop.is_set():
                self._shutdown_work()
        finally:
            self._stop.set()
            get_runtime().set_listener(None)
            self._close_listener()
            with self._conns_lock:
                connections = list(self._connections)
            for conn in connections:
                conn.close()
            try:
                (self.dir / INFO).unlink()
            except OSError:
                pass
            self._lockfile.close()          # daemon.lock first, then start.lock
            if start_lock is not None:
                start_lock.close()
        return 0

    def stop(self):
        """End :meth:`serve` from another thread, cancelling the running workflow."""
        self._stop.set()
        self._wake.set()

    def _on_signal(self, signum, _frame):
        log.info('signal %s, stopping', signal.Signals(signum).name)
        self.stop()

    def _shutdown_work(self):
        with self._lock:
            active = self._active
        if active is not None and active.future is not None:
            active.token.cancel()
            active.future.result()

    def _spawn(self, target, name, *args):
        thread = threading.Thread(target=target, args=args, name=name, daemon=True)
        thread.start()
        return thread

    def _write_info(self):
        submission._atomic(self.dir / INFO, {'pid': os.getpid(), 'host': socket.gethostname(),
                                            'started': time.time(), 'budget': self.budget})

    def _listen(self):
        path = self.dir / SOCKET
        try:
            path.unlink()                   # a stale socket: this process holds daemon.lock
        except FileNotFoundError:
            pass
        listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        listener.bind(str(path))
        listener.listen(64)
        self._listener = listener
        self._accepting = True
        self._accept_thread = self._spawn(self._accept_loop, 'tumult-daemon-accept', listener)

    def _close_listener(self):
        listener, self._listener = self._listener, None
        self._accepting = False
        if listener is None:
            return
        self._accept_thread.join()
        listener.close()
        try:
            (self.dir / SOCKET).unlink()
        except OSError:
            pass

    def _accept_loop(self, listener):
        while self._accepting:
            ready, _, _ = select.select([listener], [], [], 0.2)
            if not ready or not self._accepting:
                continue
            try:
                sock, _ = listener.accept()
            except OSError:
                continue
            conn = _Connection(sock)
            with self._conns_lock:
                self._connections.append(conn)
            self._spawn(conn.write_loop, 'tumult-daemon-write')
            self._spawn(self._read_loop, 'tumult-daemon-read', conn)

    def _idle(self):
        """``0`` when the exit rule holds, else how long to wait before checking again
        (``CHECK_S`` while there is work or a server)."""
        with self._lock:
            if self._active is not None or self._queue:
                return CHECK_S
        now = time.monotonic()
        with self._conns_lock:
            if any(conn.role == 'server' for conn in self._connections):
                return CHECK_S
            last = max([conn.accepted for conn in self._connections] + [self._last_request])
        if not self._ever_connected and now - self._started < self.grace_s:
            return self._started + self.grace_s - now
        # a client may be between connecting and submitting; it would retry on a new
        # daemon, but waiting a moment spares the restart
        return max(0.0, last + RECENT_S - now)

    def _close_if_idle(self):
        """Close the socket under ``start.lock`` if the daemon is still idle; returns the
        held lock file to release after ``daemon.lock``, or ``None`` when work or a server
        arrived meanwhile (the socket is then open again)."""
        import fcntl
        handle = open(self.dir / START_LOCK, 'a+')
        fcntl.flock(handle, fcntl.LOCK_EX)
        if self._idle():
            handle.close()
            return None
        self._close_listener()
        if self._idle():
            self._listen()
            handle.close()
            return None
        return handle

    # --- connections ------------------------------------------------------------------
    def _read_loop(self, conn):
        try:
            with conn.sock.makefile('rb') as lines:
                for line in lines:
                    if not line.strip():
                        continue
                    self._last_request = time.monotonic()
                    message = {}
                    try:
                        message = json.loads(line)
                        reply = self._handle(conn, message)
                    except Exception as exc:
                        log.warning('request %r failed', line[:200], exc_info=True)
                        reply = {'error': f'{type(exc).__name__}: {exc}'}
                    conn.send({'id': message.get('id') if isinstance(message, dict) else None,
                               **reply})
        except OSError:
            pass
        finally:
            conn.close()
            conn.sock.close()
            with self._conns_lock:
                if conn in self._connections:
                    self._connections.remove(conn)
            if conn.role == 'server':
                log.info('server pid %s disconnected', conn.pid)
            self._wake.set()

    def _handle(self, conn, message):
        op = message.get('op')
        if op == 'hello':
            return self._hello(conn, message)
        if conn.role is None:
            raise ValueError('the first request must be hello')
        if op == 'submit':
            return self._submit(Path(message['folder']))
        if op == 'cancel':
            return self._cancel(str(message['job']), message.get('folder'))
        if op == 'list':
            for root in message.get('roots') or ():
                self._register(root)
            return self._list(int(message.get('limit') or LIST_LIMIT))
        if op == 'subscribe':
            jobs = message.get('jobs')
            if jobs == 'all':
                conn.all = True
            else:
                conn.subscribed.update(str(job) for job in jobs or ())
            return {}
        if op == 'configure':
            return {'budget': self._configure(message.get('max_threads'))}
        raise ValueError(f'unknown request {op!r}')

    def _hello(self, conn, message):
        role = message.get('role')
        if role not in ROLES:
            raise ValueError(f'role must be one of {ROLES}, not {role!r}')
        conn.role, conn.pid = role, message.get('pid')
        self._ever_connected = True
        for root in message.get('roots') or ():
            self._register(root)
        if role == 'server' and message.get('max_threads') is not None:
            self._configure(message['max_threads'])
        log.info('%s pid %s connected', role, conn.pid)
        self._wake.set()        # the startup grace no longer applies
        return {'daemon_pid': os.getpid(), 'host': socket.gethostname(), 'budget': self.budget}

    def _configure(self, max_threads):
        self.budget = resolve_cpu_budget(max_threads)
        get_runtime().configure(self.budget)
        return self.budget

    def _emit(self, message, job=None):
        with self._conns_lock:
            connections = list(self._connections)
        for conn in connections:
            if conn.role is not None and (job is None or conn.wants(job)):
                conn.send(message)

    def _emit_state(self, workflow):
        summary = _summary(workflow.folder, workflow.record, workflow.threads)
        workflow.emitted = (summary['state'], summary['phase'])
        self._emit({'event': 'state', **summary}, workflow.id)
        if summary['state'] in TERMINAL:
            with self._conns_lock:
                for conn in self._connections:
                    conn.subscribed.discard(workflow.id)

    def _on_allocation(self, changed, used, budget):
        # runtime lock held: only queue messages here
        with self._lock:
            workflows = {w.id: w for w in self._queue}
            if self._active is not None:
                workflows[self._active.id] = self._active
        for owner, threads in changed.items():
            if isinstance(owner, str):
                if owner in workflows:
                    workflows[owner].threads = threads
                self._emit({'event': 'threads', 'job': owner, 'threads': threads}, owner)
        self._emit({'event': 'usage', 'used': used, 'budget': budget})

    # --- roots and records ------------------------------------------------------------
    def _stored_roots(self):
        try:
            return [Path(root) for root in json.loads((self.dir / ROOTS).read_text())]
        except (OSError, ValueError):
            return []

    def _register(self, root):
        """Know ``root`` from now on: re-queue its queued daemon records and mark running
        ones this daemon does not run interrupted."""
        root = Path(root).expanduser().resolve()
        with self._lock:
            new = root not in self._roots
            self._roots.add(root)
            if new:
                submission._atomic(self.dir / ROOTS, sorted(str(r) for r in self._roots))
        if not new:
            return
        host = socket.gethostname()
        for path in sorted((root / '_workflows').glob(f'*/{submission.RECORD}')):
            record = submission._read(path.parent)
            driver = (record or {}).get('driver') or {}
            if driver.get('kind') != 'daemon' or driver.get('host') not in (None, host):
                continue
            if record.get('state') == 'queued':
                self._submit(path.parent)
            elif record.get('state') == 'running':
                with self._lock:
                    mine = self._active is not None and self._active.folder == path.parent
                if not mine:
                    workflow = _Workflow(path.parent, record)
                    workflow.record = submission.update_record(
                        path.parent, record, state='interrupted', phase='interrupted',
                        ended=time.time(), error='the daemon running it exited')
                    log.info('workflow %s: interrupted', workflow.id)
                    self._emit_state(workflow)

    def _submit(self, folder):
        folder = Path(folder).expanduser().resolve()
        record = submission._read(folder)
        if record is None:
            raise FileNotFoundError(f'No workflow record in {folder}')
        identity = record.get('id') or folder.name
        with self._lock:
            known = (self._active is not None and self._active.id == identity) or any(
                w.id == identity for w in self._queue)
        if known:
            return {'accepted': True, 'state': record.get('state')}
        if record.get('state') != 'queued':
            return {'accepted': False, 'state': record.get('state')}
        self._register(folder.parent.parent)
        workflow = _Workflow(folder, record)
        with self._lock:
            if any(w.id == identity for w in self._queue) or (
                    self._active is not None and self._active.id == identity):
                return {'accepted': True, 'state': 'queued'}
            self._queue.append(workflow)
            self._queue.sort(key=lambda w: (w.submitted, w.id))
        self._emit_state(workflow)
        self._pump()
        return {'accepted': True, 'state': 'queued'}

    def _cancel(self, job, folder=None):
        with self._lock:
            queued = next((w for w in self._queue if w.id == job), None)
            if queued is not None:
                self._queue.remove(queued)
            active = self._active if self._active is not None and self._active.id == job else None
        if queued is not None:
            queued.record = submission.update_record(queued.folder, queued.record,
                                                     state='cancelled', phase='cancelled',
                                                     ended=time.time())
            self._emit_state(queued)
            self._wake.set()
            return {'state': 'cancelled'}
        if active is not None:
            active.token.cancel()
            return {'state': 'cancelling'}
        external = self._external.get(job)
        if external is not None:
            if jobrecord.cancel(external['folder']):
                self._external.pop(job, None)
                return {'state': 'cancelled'}
            return {'state': jobrecord.state_of(external['folder'])}
        record = submission._read(folder) if folder else None
        if record is None:
            raise KeyError(f'unknown job {job!r}')
        driver = record.get('driver') or {}
        if record.get('state') == 'queued' and driver.get('kind') == 'daemon':
            # never submitted to this daemon: settle it here rather than leave it queued
            record = submission.update_record(folder, record, state='cancelled',
                                              phase='cancelled', ended=time.time())
            self._emit_state(_Workflow(folder, record))
        return {'state': record.get('state')}

    def _list(self, limit):
        with self._lock:
            live = ([self._active] if self._active is not None else []) + list(self._queue)
            roots = sorted(self._roots)
        jobs = [_summary(w.folder, w.record, w.threads) for w in live]
        seen = {w.id for w in live}
        for root in roots:
            settled = []
            for path in (root / '_workflows').glob(f'*/{submission.RECORD}'):
                record = submission._read(path.parent)
                if record is None or (record.get('id') or path.parent.name) in seen:
                    continue
                if (record.get('driver') or {}).get('kind') == 'slurm' or \
                        record.get('state') not in TERMINAL:
                    record['state'] = submission.WorkflowHandle(path.parent).status
                settled.append(_summary(path.parent, record))
            settled.sort(key=lambda job: job.get('submitted') or 0, reverse=True)
            jobs.extend(settled[:limit])
        jobs.extend(self._external.values())
        runtime = get_runtime()
        return {'jobs': jobs, 'budget': runtime.budget, 'used': runtime.used}

    # --- running workflows -------------------------------------------------------------
    def _pump(self):
        with self._lock:
            if self._active is not None or not self._queue or self._stop.is_set():
                return
            workflow = self._queue.pop(0)
            self._active = workflow
        try:
            self._start(workflow)
        except BaseException:
            with self._lock:
                self._active = None
            raise

    def _start(self, workflow):
        try:
            kind, args = submission.decode_task(workflow.record)
        except Exception as exc:
            log.exception('workflow %s has an unreadable request', workflow.id)
            workflow.record = submission.update_record(
                workflow.folder, workflow.record, state='failed', phase='failed',
                ended=time.time(), error=f'unreadable request: {exc}')
            self._settled(workflow)
            return

        def claim(current):
            driver = current.get('driver') or {}
            driver.update(kind='daemon', pid=os.getpid(), host=socket.gethostname(),
                          started=time.time(), budget=self.budget)
            current['driver'] = driver

        workflow.record = submission.mutate_record(workflow.folder, claim)
        workflow.token = CancelToken()

        def changed(record, progress_stamp=None):
            if progress_stamp is not None:
                self._emit({'event': 'progress', 'job': workflow.id, 'stamp': progress_stamp},
                           workflow.id)
                return
            workflow.record = record
            if workflow.emitted != (record.get('state'), record.get('phase')):
                self._emit_state(workflow)

        def work():
            try:
                submission.execute_recorded(workflow.folder, kind, args, workflow.token,
                                            workflow.record, listener=changed)
            except (CancelledError, Exception):
                pass        # the record names the outcome; the traceback is in the log
            finally:
                self._settled(workflow)

        log.info('workflow %s: started with budget %d', workflow.id, self.budget)
        workflow.future = WORKFLOWS.submit(work)

    def _settled(self, workflow):
        log.info('workflow %s: %s', workflow.id, workflow.record.get('state'))
        with self._lock:
            if self._active is workflow:
                self._active = None
        record = workflow.record
        if workflow.emitted != (record.get('state'), record.get('phase')):
            self._emit_state(workflow)
        self._pump()
        self._wake.set()

    # --- external runs -------------------------------------------------------------------
    def _scan_loop(self):
        while not self._stop.wait(SCAN_S):
            with self._conns_lock:
                serving = any(conn.role == 'server' for conn in self._connections)
            if serving:
                try:
                    self.scan_external()
                except Exception:
                    log.exception('external-run scan failed')

    def scan_external(self):
        """Refresh the live dispatch runs under the registered roots that no workflow
        record names; emits ``external`` when the set changed."""
        with self._lock:
            roots = sorted(self._roots)
        found = {}
        cache = {}
        for root in roots:
            accounted = set()
            for path in (root / '_workflows').glob(f'*/{submission.RECORD}'):
                accounted.update(str(Path(run)) for run in
                                 ((submission._read(path.parent) or {}).get('runs') or ()))
            folders = {path.parent for name in (jobrecord.JOB_NAME, jobrecord.NATIVE_NAME)
                       for path in root.rglob(name) if '_workflows' not in path.parts}
            for folder in folders:
                if str(folder) in accounted:
                    continue
                stamp = max((p.stat().st_mtime_ns for p in
                             (folder / jobrecord.JOB_NAME, folder / jobrecord.NATIVE_NAME)
                             if p.exists()), default=0)
                cached = self._external_cache.get(folder)
                record = cached[1] if cached and cached[0] == stamp else jobrecord.read(folder)
                cache[folder] = (stamp, record)
                if record is None or record.get('state') not in ('running', 'submitted'):
                    continue            # a settled record: nothing to re-resolve
                record = dict(record, resolved=jobrecord.state_of(folder, record))
                if record['resolved'] in ('running', 'submitted'):
                    summary = _external_summary(record)
                    found[summary['job']] = summary
        self._external_cache = cache
        if found != self._external:
            self._external = found
            self._emit({'event': 'external', 'jobs': list(found.values())})
        return found


# --- the client -----------------------------------------------------------------------------
class DaemonError(RuntimeError):
    """The daemon refused a request; the message is its error."""


def _connect_socket(path):
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        sock.connect(str(path))
    except OSError:
        sock.close()
        return None
    return sock


def _start_daemon(base):
    """Spawn ``python -m tumult.synthesis.daemon`` in its own session, logging to
    ``daemon.log``. Relative paths in requests resolve against the working directory it
    inherits from this process."""
    env = {**os.environ}
    package_root = str(Path(__file__).resolve().parents[2])
    env['PYTHONPATH'] = os.pathsep.join(
        part for part in (package_root, env.get('PYTHONPATH')) if part)
    with open(base / LOGFILE, 'ab') as log_file:
        subprocess.Popen([sys.executable, '-m', 'tumult.synthesis.daemon'], cwd=os.getcwd(),
                         env=env, stdin=subprocess.DEVNULL, stdout=log_file,
                         stderr=subprocess.STDOUT, start_new_session=True)


def open_socket(start=True):
    """A connected socket to the daemon, starting one when none listens and ``start``."""
    if os.name != 'posix':
        raise OSError('The workflow daemon needs a POSIX host')
    base = directory()
    path = base / SOCKET
    sock = _connect_socket(path)
    if sock is not None or not start:
        if sock is None:
            raise ConnectionError(f'No workflow daemon listens on {path}')
        return sock
    with _flocked(base / START_LOCK):
        sock = _connect_socket(path)
        if sock is None:
            _start_daemon(base)
            deadline = time.monotonic() + START_TIMEOUT_S
            while sock is None and time.monotonic() < deadline:
                time.sleep(0.02)
                sock = _connect_socket(path)
    if sock is None:
        raise ConnectionError(f'The workflow daemon did not start; see {base / LOGFILE}')
    return sock


class Client:
    """A connection to the daemon with request/reply, subscriptions and event dispatch.

    Events are delivered on one dispatcher thread, never the reader, so a callback may
    make requests of its own. ``role='server'`` keeps the daemon alive while connected
    and reconnects (starting a daemon) as soon as the connection drops; a client role
    reconnects on its next request. Either way the hello and the subscriptions are sent
    again, and every callback receives ``{'event': 'reconnect'}``, since events may have
    been missed meanwhile. Requests are idempotent, so one interrupted by a dropped
    connection is sent once more on the new one."""

    def __init__(self, role='client', *, roots=(), max_threads=None, start=True):
        if role not in ROLES:
            raise ValueError(f'role must be one of {ROLES}, not {role!r}')
        self.role = role
        self.roots = [str(Path(root).expanduser().resolve()) for root in roots]
        self.max_threads = max_threads
        self.start = start
        self.info = {}
        self._ids = itertools.count(1)
        self._lock = threading.RLock()          # short critical sections only
        self._connect_lock = threading.Lock()   # one connection attempt at a time
        self._send_lock = threading.Lock()
        self._pending = {}
        self._followers = {}        # job id -> [callback]
        self._listeners = []        # callbacks for every event
        self._all = False
        self._sock = None
        self._closed = False
        self._events = queue.Queue()
        threading.Thread(target=self._dispatch_loop, name='tumult-daemon-events',
                         daemon=True).start()
        self._connect()

    # --- connection ------------------------------------------------------------------
    @property
    def connected(self):
        return self._sock is not None

    def _connect(self):
        """Connect, say hello and subscribe again unless connected; returns the socket."""
        with self._connect_lock:
            sock = self._sock
            if sock is not None:
                return sock
            if self._closed:
                raise ConnectionError('The daemon client is closed')
            sock = open_socket(self.start)
            threading.Thread(target=self._read_loop, args=(sock,), name='tumult-daemon-reply',
                             daemon=True).start()
            reconnect = bool(self.info)
            try:
                info = self._request_on(sock, 'hello', timeout=START_TIMEOUT_S, role=self.role,
                                        pid=os.getpid(), host=socket.gethostname(),
                                        roots=self.roots, max_threads=self.max_threads)
                with self._lock:
                    jobs = 'all' if self._all else sorted(self._followers)
                if jobs:
                    self._request_on(sock, 'subscribe', timeout=START_TIMEOUT_S, jobs=jobs)
            except BaseException:
                sock.close()
                raise
            self.info = info
            with self._lock:
                self._sock = sock
        if reconnect:
            self._events.put({'event': 'reconnect'})
        return sock

    def _read_loop(self, sock):
        try:
            with sock.makefile('rb') as lines:
                for line in lines:
                    message = json.loads(line)
                    if 'event' in message:
                        self._events.put(message)
                        continue
                    with self._lock:
                        waiter = self._pending.pop(message.get('id'), None)
                    if waiter is not None:
                        waiter[1] = message
                        waiter[0].set()
        except (OSError, ValueError):
            pass
        self._dropped(sock)

    def _dropped(self, sock):
        with self._lock:
            current = self._sock is sock
            if current:
                self._sock = None
            pending = [waiter for waiter in self._pending.values() if waiter[2] is sock]
            self._pending = {key: waiter for key, waiter in self._pending.items()
                             if waiter[2] is not sock}
        try:
            sock.close()
        except OSError:
            pass
        for waiter in pending:
            waiter[0].set()             # no reply: the requester sends again
        if current and self.role == 'server' and not self._closed:
            threading.Thread(target=self._reconnect_loop, name='tumult-daemon-reconnect',
                             daemon=True).start()

    def _reconnect_loop(self):
        delay = 0.1
        while not self._closed and self._sock is None:
            try:
                self._connect()
                log.info('reconnected to the workflow daemon')
                return
            except (OSError, DaemonError):
                log.warning('workflow daemon unreachable; retrying in %.1f s', delay,
                            exc_info=True)
                time.sleep(delay)
                delay = min(5.0, delay * 2)

    def close(self):
        """Close the connection; a server role stops keeping the daemon alive."""
        with self._lock:
            self._closed = True
            sock, self._sock = self._sock, None
        if sock is not None:
            try:
                sock.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass
            sock.close()
        self._events.put(None)

    # --- requests --------------------------------------------------------------------
    def _request_on(self, sock, op, timeout=None, **args):
        identity = next(self._ids)
        waiter = [threading.Event(), None, sock]
        with self._lock:
            self._pending[identity] = waiter
        data = (json.dumps({'id': identity, 'op': op, **args}, default=str) + '\n').encode()
        try:
            with self._send_lock:
                sock.sendall(data)
        except OSError as exc:
            with self._lock:
                self._pending.pop(identity, None)
            raise ConnectionError(f'workflow daemon connection lost: {exc}') from exc
        if not waiter[0].wait(timeout):
            with self._lock:
                self._pending.pop(identity, None)
            raise TimeoutError(f'the workflow daemon did not answer {op}')
        reply = waiter[1]
        if reply is None:
            raise ConnectionError(f'workflow daemon connection lost during {op}')
        if 'error' in reply:
            raise DaemonError(reply['error'])
        return reply

    def request(self, op, *, timeout=60, **args):
        """Send one request and return its reply; raises :class:`DaemonError` for a
        refused one and ``ConnectionError`` when no daemon can be reached."""
        for attempt in (0, 1):
            sock = self._connect()
            try:
                return self._request_on(sock, op, timeout=timeout, **args)
            except ConnectionError:
                self._dropped(sock)
                if attempt or self._closed:
                    raise

    # --- events ------------------------------------------------------------------------
    def listen(self, callback):
        """Deliver every event to ``callback(event)`` and subscribe to all jobs."""
        with self._lock:
            self._listeners.append(callback)
            self._all = True
        self.request('subscribe', jobs='all')

    def follow(self, job, callback):
        """Deliver the events of ``job`` to ``callback(event)``; returns a function that
        stops following."""
        with self._lock:
            self._followers.setdefault(job, []).append(callback)
        self.request('subscribe', jobs=[job])

        def unfollow():
            with self._lock:
                callbacks = self._followers.get(job, [])
                if callback in callbacks:
                    callbacks.remove(callback)
                if not callbacks:
                    self._followers.pop(job, None)
        return unfollow

    def _dispatch_loop(self):
        while True:
            event = self._events.get()
            if event is None:
                return
            with self._lock:
                job = event.get('job')
                if event['event'] == 'reconnect':
                    callbacks = [cb for cbs in self._followers.values() for cb in cbs]
                else:
                    callbacks = list(self._followers.get(job, ())) if job else []
                callbacks += self._listeners
            for callback in callbacks:
                try:
                    callback(event)
                except Exception:
                    log.exception('workflow daemon event callback failed for %s', event)


LIVE = ('queued', 'running')


class LiveJobs:
    """The daemon's running and queued workflows, its external runs and the budget in
    use, kept current from the events of ``client``. Watchers are called with each event
    after it has been applied, on the client's event thread."""

    def __init__(self, client):
        self.client = client
        self.used, self.budget = 0, None
        self._jobs = {}
        self._lock = threading.Lock()
        self._watchers = []
        client.listen(self._on_event)
        self.refresh()

    def refresh(self, roots=(), limit=LIST_LIMIT):
        """Replace the view with a fresh ``list`` (registering ``roots``); returns the
        reply, whose settled records the view itself does not keep."""
        reply = self.client.request('list', roots=[str(root) for root in roots], limit=limit)
        with self._lock:
            self._jobs = {job['job']: job for job in reply['jobs']
                          if job['type'] == 'external' or job['state'] in LIVE}
            self.used, self.budget = reply['used'], reply['budget']
        self._notify({'event': 'refresh'})
        return reply

    def jobs(self):
        with self._lock:
            return [dict(job) for job in self._jobs.values()]

    def get(self, job):
        with self._lock:
            found = self._jobs.get(job)
            return dict(found) if found is not None else None

    def external_threads(self):
        """Threads of the live external runs on this host, outside the budget."""
        with self._lock:
            return sum(job.get('threads') or 0 for job in self._jobs.values()
                       if job['type'] == 'external')

    def watch(self, callback):
        """Call ``callback(event)`` after every change; returns a function that stops."""
        with self._lock:
            self._watchers.append(callback)

        def unwatch():
            with self._lock:
                if callback in self._watchers:
                    self._watchers.remove(callback)
        return unwatch

    def _on_event(self, event):
        kind = event['event']
        if kind == 'reconnect':
            self.refresh()
            return
        with self._lock:
            if kind == 'state':
                job = {key: value for key, value in event.items() if key != 'event'}
                if job['state'] in LIVE:
                    self._jobs[job['job']] = {**self._jobs.get(job['job'], {}), **job}
                else:
                    self._jobs.pop(job['job'], None)
            elif kind == 'threads' and event['job'] in self._jobs:
                self._jobs[event['job']]['threads'] = event['threads']
            elif kind == 'usage':
                self.used, self.budget = event['used'], event['budget']
            elif kind == 'external':
                self._jobs = {key: job for key, job in self._jobs.items()
                              if job['type'] != 'external'}
                self._jobs.update((job['job'], job) for job in event['jobs'])
        self._notify(event)

    def _notify(self, event):
        with self._lock:
            watchers = list(self._watchers)
        for callback in watchers:
            try:
                callback(event)
            except Exception:
                log.exception('live job watcher failed')


_shared = None
_live = None
_shared_lock = threading.Lock()


def client():
    """This process's connection to the daemon, connecting as ``client`` on first use.
    A process that called :func:`connect_server` gets its server connection instead."""
    global _shared
    with _shared_lock:
        if _shared is None or _shared._closed:
            _shared = Client('client')
        return _shared


def current():
    """This process's connection if it has one, without connecting."""
    with _shared_lock:
        return _shared if _shared is not None and not _shared._closed else None


def live_jobs():
    """The :class:`LiveJobs` view of this process's connection, or ``None`` while the
    process has none (a document built outside a Bokeh server, say)."""
    global _live
    connection = current()
    if connection is None:
        return None
    with _shared_lock:
        if _live is None or _live.client is not connection:
            _live = None
    if _live is None:
        view = LiveJobs(connection)
        with _shared_lock:
            if _live is None or _live.client is not connection:
                _live = view
    return _live


def connect_server(roots=(), max_threads=None):
    """Make this process's connection a ``server`` one that keeps the daemon alive and
    registers ``roots``; returns it. Used once per Bokeh server process."""
    global _shared
    with _shared_lock:
        previous, _shared = _shared, Client('server', roots=roots, max_threads=max_threads)
    if previous is not None:
        previous.close()
    return _shared


def disconnect():
    """Close this process's connection, if any."""
    global _shared, _live
    with _shared_lock:
        previous, _shared, _live = _shared, None, None
    if previous is not None:
        previous.close()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) > 1:
        print(__doc__, file=sys.stderr)
        return 3
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    return Daemon(int(argv[0]) if argv else 0).serve()


if __name__ == '__main__':
    sys.exit(main())
