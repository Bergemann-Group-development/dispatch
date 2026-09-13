"""One long-lived worker per workflow root that runs queued workflows off the server.

    python -m tumult.synthesis.daemon <output root>

The daemon serves ``<root>/_workflows``: every second it claims records that are
``queued`` for driver kind ``daemon`` and runs them on its own ``WORKFLOWS`` pool with
the shared :class:`runtime.ExecutionService`, so concurrent workflows share identical
dispatch runs and one thread budget as they do inside the server. It cancels a running
workflow when its record carries ``cancel_requested`` (written by any handle), ends
everything on SIGTERM/SIGINT/SIGHUP, and exits on its own after ``IDLE_S`` without work.
The submitter starts it when none is alive; ``daemon.json`` beside the records names the
pid and host, and ``daemon.lock`` is held for its lifetime so a second start finds the
lock taken and leaves.

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
from .request import ResolvedSynthRequest
from .runtime import WORKFLOWS
from . import submission

log = logging.getLogger(__name__)

PIDFILE = 'daemon.json'
LOCKFILE = 'daemon.lock'
LOGFILE = 'daemon.log'
POLL_S = 1.0
IDLE_S = 600.0


class Daemon:
    def __init__(self, root, idle_s=IDLE_S):
        self.base = Path(root).expanduser().resolve() / '_workflows'
        self.idle_s = idle_s
        self._tokens = {}           # workflow id -> CancelToken of a running workflow
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
            'root': str(self.base.parent)})
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
        log.info('serving %s as pid %d on %s', self.base, os.getpid(), socket.gethostname())
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
            elif record.get('state') == 'queued' and not (record.get('driver') or {}).get('pid'):
                self._claim(folder, record)
        with self._lock:
            for identity, future in list(self._futures.items()):
                if future.done():
                    self._futures.pop(identity)
                    self._tokens.pop(identity)
                    self._last_work = time.monotonic()

    def _claim(self, folder, record):
        identity = record.get('id') or folder.name
        try:
            resolved = ResolvedSynthRequest.from_dict(record['request'])
        except (KeyError, ValueError, TypeError) as exc:
            submission.update_record(folder, record, state='failed', phase='failed',
                                     ended=time.time(), error=f'unreadable request: {exc}')
            log.exception('workflow %s has an unreadable request', identity)
            return

        def claim(current):
            driver = current.get('driver') or {}
            if driver.get('pid') or current.get('state') != 'queued':
                raise _Taken()
            driver.update(pid=os.getpid(), host=socket.gethostname(), started=time.time())
            current['driver'] = driver

        try:
            submission.mutate_record(folder, claim)
        except _Taken:
            return
        if record.get('cancel_requested'):
            submission.update_record(folder, state='cancelled', phase='cancelled',
                                     ended=time.time())
            return
        token = CancelToken()

        def work():
            try:
                submission.execute_recorded(folder, resolved, token, record)
            except (CancelledError, Exception):
                pass        # the record names the outcome; the traceback is in the log

        with self._lock:
            self._tokens[identity] = token
            self._futures[identity] = WORKFLOWS.submit(work)
            self._last_work = time.monotonic()
        log.info('workflow %s: claimed', identity)


class _Taken(Exception):
    """Another daemon claimed the record between the read and the locked write."""


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


def ensure(root):
    """Start a daemon for ``root`` unless one is alive; returns its pid or ``None`` when
    the start was only requested (the new daemon writes its pidfile once it holds the
    lock). Never blocks on the child."""
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
        subprocess.Popen([sys.executable, '-m', 'tumult.synthesis.daemon', str(base.parent)],
                         cwd=os.getcwd(), env=env, stdin=subprocess.DEVNULL,
                         stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True)
    return None


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 3
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    return Daemon(argv[0]).serve()


if __name__ == '__main__':
    sys.exit(main())
