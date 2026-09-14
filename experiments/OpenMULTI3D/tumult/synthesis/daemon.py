"""One long-lived worker per workflow root that runs queued workflows off the server.

    python -m tumult.synthesis.daemon <output root> [max_threads]

The daemon serves ``<root>/_workflows``: every second it claims and runs at most one
queued record for driver kind ``daemon`` at a time (FIFO by submission time, workflow id
as tie-breaker), leaving every other queued record untouched until the active one
settles (see "Serialize workflows in the daemon" in ``plan/multi-processing.md``). This
gives the active workflow the daemon's full, once-resolved CPU budget for every
parallelizable native phase instead of splitting it between concurrent workflows. It
cancels the active workflow when its record carries ``cancel_requested`` (written by any
handle), settles a *queued* cancellation without ever claiming it, ends everything on
SIGTERM/SIGINT/SIGHUP, and exits on its own after ``IDLE_S`` without work. The submitter
starts it when none is alive; ``daemon.json`` beside the records names the pid, host and
resolved budget, and ``daemon.lock`` is held for its lifetime so a second start finds the
lock taken and leaves.

``max_threads`` (default ``0``) is resolved once at startup with
:func:`runtime.resolve_cpu_budget`: ``0`` reserves one logical CPU
(``max(1, cpu_count - 1)``); an explicit value opts into that reserved CPU but is capped
at the detected logical CPU count. The resolved budget configures the shared
:class:`runtime.ExecutionService` once and is copied into every claimed workflow's
``driver`` block for observability.

Interpreter startup is paid once per daemon lifetime rather than once per workflow, which
is the point of this over ``tumult.synthesis.driver``: short interactive syntheses run
detached without a visible delay. A workflow that was queued when its daemon died is
picked up by the next one; a workflow that was running is reported ``interrupted``.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
import signal
import socket
import sys
import threading
import time

from .cancel import CancelToken, CancelledError
from .runtime import WORKFLOWS, get_runtime, resolve_cpu_budget
from . import submission

log = logging.getLogger(__name__)

PIDFILE = 'daemon.json'
LOCKFILE = 'daemon.lock'
LOGFILE = 'daemon.log'
POLL_S = 1.0
IDLE_S = 600.0


class Daemon:
    def __init__(self, root, idle_s=IDLE_S, max_threads=0):
        self.base = Path(root).expanduser().resolve() / '_workflows'
        self.idle_s = idle_s
        self.budget = resolve_cpu_budget(max_threads)
        self._tokens = {}           # workflow id -> CancelToken of the active workflow
        self._futures = {}
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._last_work = time.monotonic()

    # --- lifetime ------------------------------------------------------------------
    def _acquire(self):
        import fcntl
        self.base.mkdir(parents=True, exist_ok=True)
        self._lockfile = open(self.base / LOCKFILE, 'a+')
        try:
            fcntl.flock(self._lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return False
        submission._atomic(self.base / PIDFILE, {
            'pid': os.getpid(), 'host': socket.gethostname(), 'started': time.time(),
            'root': str(self.base.parent), 'budget': self.budget})
        return True

    def _release(self):
        try:
            (self.base / PIDFILE).unlink()
        except OSError:
            pass
        self._lockfile.close()

    def serve(self):
        if not self._acquire():
            log.info('another daemon holds %s', self.base / LOCKFILE)
            return 0
        if threading.current_thread() is threading.main_thread():
            for name in ('SIGTERM', 'SIGINT', 'SIGHUP'):
                if hasattr(signal, name):
                    signal.signal(getattr(signal, name), self._on_signal)
        # Configure the shared execution service once, from this daemon's resolved
        # budget; the active workflow must not mutate it while work is running
        # (plan/multi-processing.md, "CPU budget").
        get_runtime().configure(self.budget)
        log.info('serving %s as pid %d on %s with budget %d',
                self.base, os.getpid(), socket.gethostname(), self.budget)
        try:
            while not self._stop.is_set():
                self._tick()
                if not self._futures and time.monotonic() - self._last_work > self.idle_s:
                    log.info('idle for %.0f s, exiting', self.idle_s)
                    break
                self._stop.wait(POLL_S)
            if self._stop.is_set():
                self._cancel_all()
            for future in list(self._futures.values()):
                try:
                    future.result()
                except BaseException:
                    pass
        finally:
            self._release()
        return 0

    def _on_signal(self, signum, _frame):
        log.info('signal %s, stopping', signal.Signals(signum).name)
        self._stop.set()

    def _cancel_all(self):
        with self._lock:
            tokens = list(self._tokens.values())
        for token in tokens:
            token.cancel()

    # --- one poll ---------------------------------------------------------------------
    def _tick(self):
        """Reap the active workflow if it settled, service cancellation for it, settle
        every other queued record that was cancelled while waiting (without ever
        claiming it), and otherwise claim at most the single oldest still-queued record
        (FIFO by submission time, workflow id as tie-breaker). A queued record nobody
        cancels is left untouched: no pid, no coordinator, so it stays cleanly
        recoverable by the next daemon."""
        self._reap()
        candidates = []
        for path in sorted(self.base.glob(f'*/{submission.RECORD}')):
            folder = path.parent
            record = submission._read(folder)
            if record is None or (record.get('driver') or {}).get('kind') != 'daemon':
                continue
            identity = record.get('id') or folder.name
            with self._lock:
                mine = identity in self._tokens
            if mine:
                if record.get('cancel_requested') and record.get('state') == 'running':
                    log.info('workflow %s: cancel requested', identity)
                    self._tokens[identity].cancel()
                continue
            if record.get('state') != 'queued' or (record.get('driver') or {}).get('pid'):
                continue
            if record.get('cancel_requested'):
                # Settle it here rather than waiting for its FIFO turn: a cancelled
                # queued workflow must never be claimed, even behind an active one.
                self._settle_queued_cancel(folder)
                continue
            candidates.append((record.get('submitted', 0), identity, folder))
        with self._lock:
            active = bool(self._tokens)
        if active or not candidates:
            return
        candidates.sort(key=lambda item: (item[0], item[1]))
        _, _, folder = candidates[0]
        record = submission._read(folder)
        if record is None:
            return
        self._claim(folder, record)

    def _settle_queued_cancel(self, folder):
        def settle(current):
            if current.get('state') != 'queued' or (current.get('driver') or {}).get('pid'):
                raise _Taken()      # claimed or settled elsewhere between the read and here
            current.update(state='cancelled', phase='cancelled', ended=time.time())
        try:
            submission.mutate_record(folder, settle)
        except _Taken:
            pass

    def _reap(self):
        with self._lock:
            for identity, future in list(self._futures.items()):
                if future.done():
                    self._futures.pop(identity)
                    self._tokens.pop(identity)
                    self._last_work = time.monotonic()

    def _claim(self, folder, record):
        identity = record.get('id') or folder.name
        try:
            kind, args = submission.decode_task(record)
        except Exception as exc:
            submission.update_record(folder, record, state='failed', phase='failed',
                                     ended=time.time(), error=f'unreadable request: {exc}')
            log.exception('workflow %s has an unreadable request', identity)
            return

        def claim(current):
            driver = current.get('driver') or {}
            if driver.get('pid') or current.get('state') != 'queued':
                raise _Taken()
            if current.get('cancel_requested'):
                # Settle atomically with the claim attempt: a cancelled-while-queued
                # record must never get a pid (plan: "Cancelling an unclaimed workflow
                # ... settles it ... without launching native work").
                current.update(state='cancelled', phase='cancelled', ended=time.time())
                raise _Cancelled()
            driver.update(pid=os.getpid(), host=socket.gethostname(), started=time.time(),
                          budget=self.budget)
            current['driver'] = driver

        try:
            submission.mutate_record(folder, claim)
        except (_Taken, _Cancelled):
            return
        token = CancelToken()

        def work():
            try:
                submission.execute_recorded(folder, kind, args, token, record)
            except (CancelledError, Exception):
                pass        # the record names the outcome; the traceback is in the log

        with self._lock:
            self._tokens[identity] = token
            self._futures[identity] = WORKFLOWS.submit(work)
            self._last_work = time.monotonic()
        log.info('workflow %s: claimed with budget %d', identity, self.budget)


class _Taken(Exception):
    """Another daemon claimed the record between the read and the locked write."""


class _Cancelled(Exception):
    """The record was cancelled while queued; ``claim`` settled it, so no worker runs."""


# --- the client side --------------------------------------------------------------------
def status(root):
    """``(pid, host)`` of the live daemon serving ``root``, else ``None``. The pid counts
    only when it still runs this module for that root on this host."""
    base = Path(root).expanduser().resolve() / '_workflows'
    try:
        import json
        info = json.loads((base / PIDFILE).read_text())
    except (OSError, ValueError):
        return None
    if info.get('host') != socket.gethostname():
        return None
    if submission._process_alive(info.get('pid'), 'tumult.synthesis.daemon') is False:
        return None
    return info.get('pid'), info.get('host')


def ensure(root, max_threads=0):
    """Start a daemon for ``root`` unless one is alive; returns its pid or ``None`` when
    the start was only requested (the new daemon writes its pidfile once it holds the
    lock). Never blocks on the child.

    ``max_threads`` seeds the CPU budget of a newly started daemon (see
    :func:`runtime.resolve_cpu_budget`); it is the submission's own Engine preference,
    passed here because it is otherwise unavailable once the daemon is already running
    from an earlier submission's choice. A live daemon keeps whatever budget it already
    resolved at its own startup."""
    import subprocess
    if os.name != 'posix':
        raise ValueError('The workflow daemon needs a POSIX host')
    live = status(root)
    if live is not None:
        return live[0]
    base = Path(root).expanduser().resolve() / '_workflows'
    base.mkdir(parents=True, exist_ok=True)
    env = {**os.environ}
    package_root = str(Path(__file__).resolve().parents[2])
    env['PYTHONPATH'] = os.pathsep.join(
        part for part in (package_root, env.get('PYTHONPATH')) if part)
    with open(base / LOGFILE, 'ab') as log_file:
        subprocess.Popen([sys.executable, '-m', 'tumult.synthesis.daemon', str(base.parent),
                          str(int(max_threads or 0))],
                         cwd=os.getcwd(), env=env, stdin=subprocess.DEVNULL,
                         stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True)
    return None


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) not in (1, 2):
        print(__doc__, file=sys.stderr)
        return 3
    max_threads = int(argv[1]) if len(argv) == 2 else 0
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    return Daemon(argv[0], max_threads=max_threads).serve()


if __name__ == '__main__':
    sys.exit(main())
