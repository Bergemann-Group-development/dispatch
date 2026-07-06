"""Cooperative cancellation for in-flight syntheses, with hard process-group kill.

Extracted from the synthesis monolith in the Phase 6 decomposition: everything about
stopping a run lives here — the token the GUI holds, the exception it raises, and the
process-group plumbing that lets a cancel kill ``dispatch.x`` children too.
"""
from __future__ import annotations

import os
import signal
import threading


class CancelledError(Exception):
    """Raised inside a synthesis when its :class:`CancelToken` has been cancelled.

    Deliberately *not* a subclass of ``ValueError``/``RuntimeError`` so the fit's
    per-point ``except (ValueError, RuntimeError)`` guard does not swallow it — a cancel
    must abort the whole run, not just drop one grid point.
    """


def _init_pool_worker():
    """ProcessPoolExecutor initializer: put each worker (and the ``dispatch.x`` child it
    later launches via ``subprocess.run``) into its own process group, so cancelling a run
    can hard-kill the whole subtree with one ``killpg`` — killing only the worker would
    orphan a still-running ``dispatch.x`` (reparented to init), which is what we must avoid.
    """
    if hasattr(os, "setsid"):
        try:
            os.setsid()
        except OSError:
            pass


def _kill_process_group(pid):
    """SIGKILL the process group led by ``pid`` (worker + its ``dispatch.x`` child).

    Guards against the brief window after a worker spawns but before its initializer has
    run ``setsid``: while its pgid still equals *our* group we kill only the bare pid, so a
    cancel can never signal the main application's process group.
    """
    if pid is None or not hasattr(os, "killpg"):
        return
    try:
        pgid = os.getpgid(pid)
    except (ProcessLookupError, PermissionError):
        return
    try:
        if pgid != os.getpgid(0):       # worker established its own group: kill the tree
            os.killpg(pgid, signal.SIGKILL)
        else:                            # pre-setsid race: no descendants yet, kill the pid
            os.kill(pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


class CancelToken:
    """Cooperative cancellation for one in-flight synthesis, with hard process kill.

    The GUI holds one token per running synthesis. :meth:`cancel` flips a flag (checked at
    safe points via :meth:`check`) and, if a ``dispatch`` process pool is live, kills every
    worker process group — the workers and their ``dispatch.x`` children. The synthesis runs
    in a worker *thread* of the same process, so the pool object (and its worker PIDs) are
    reachable here; ``cancel`` is called from the Bokeh document thread, hence the lock.
    """

    def __init__(self):
        self._cancelled = False
        self._pool = None
        self._lock = threading.Lock()

    @property
    def cancelled(self):
        return self._cancelled

    def check(self):
        if self._cancelled:
            raise CancelledError()

    def set_pool(self, executor):
        """Register the live dispatch pool so a cancel can kill its workers."""
        with self._lock:
            self._pool = executor
            if self._cancelled:
                self._kill_pool_locked()   # cancelled during pool setup — kill at once

    def clear_pool(self):
        with self._lock:
            self._pool = None

    def cancel(self):
        with self._lock:
            self._cancelled = True
            self._kill_pool_locked()

    def _kill_pool_locked(self):
        pool = self._pool
        if pool is None:
            return
        # ProcessPoolExecutor keeps its live workers in ``_processes`` {pid: Process};
        # it lives in this (main) process, so the PIDs are ours to signal.
        for proc in list(getattr(pool, "_processes", {}).values()):
            _kill_process_group(proc.pid)
