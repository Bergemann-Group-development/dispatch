"""Shared, budgeted execution of dispatch jobs within one Python process.

Workflow coordinators never run on this executor. They submit bounded batches;
this service schedules their leaf jobs round-robin and shares identical in-flight
outputs. A subscriber cancelling a shared job does not cancel other subscribers.
"""
from __future__ import annotations

from collections import OrderedDict, deque
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from contextlib import contextmanager
from contextvars import ContextVar
import logging
import os
import threading
import time

from .cancel import CancelToken, CancelledError


@dataclass(eq=False)
class _Task:
    key: str
    signature: str
    owner: object
    fn: object
    threads: int
    token: CancelToken = field(default_factory=CancelToken)
    subscribers: dict = field(default_factory=dict)
    running: bool = False
    allocated: int = 0


def resolve_cpu_budget(max_threads=0):
    """The daemon/root CPU budget, resolved once from the "CPU budget" rule of
    ``plan/multi-processing.md``: worker count, treating one logical CPU as one slot.

    ``max_threads`` falsy (``0``/``None``, the default) reserves one logical CPU by
    returning ``max(1, cpu_count - 1)``; an explicit positive value opts into the
    otherwise reserved CPU but is still capped at the detected logical CPU count. A
    host that cannot report its CPU count is treated as having exactly one.
    """
    cpu = os.cpu_count()
    if not cpu or cpu < 1:
        return 1
    max_threads = int(max_threads or 0)
    if max_threads > 0:
        return max(1, min(max_threads, cpu))
    return max(1, cpu - 1)


class ExecutionService:
    """Persistent supervisors with a shared CPU budget and per-job cancellation.

    ``fn(token, threads)`` executes outside the scheduler lock. Requests with the
    same key and signature share a result; callers must treat it as read-only.
    Lowering the budget lets existing jobs drain before admitting more work.
    """

    def __init__(self, max_threads=None):
        cpu = os.cpu_count() or 1
        # An explicit budget is a ceiling, never an excuse to exceed the machine's
        # logical CPU count (plan/multi-processing.md, "CPU budget").
        self._capacity = min(int(max_threads), cpu) if max_threads else cpu
        self._capacity = max(1, self._capacity)
        self._budget = self._capacity
        self._executor = ThreadPoolExecutor(max_workers=self._capacity,
                                            thread_name_prefix="dispatch")
        self._lock = threading.RLock()
        self._queues = OrderedDict()
        self._tasks = {}
        self._running = set()
        self._used = 0
        self._closed = False

    def configure(self, max_threads):
        with self._lock:
            self._budget = min(self._capacity, max(1, int(max_threads)))
            self._pump()

    @property
    def budget(self):
        with self._lock:
            return self._budget

    def submit(self, key, signature, owner, fn, *, threads=1, token=None):
        if token is not None:
            token.check()
        future = Future()
        with self._lock:
            if self._closed:
                raise RuntimeError("Execution service is shut down")
            task = self._tasks.get(key)
            if task is not None and task.signature != signature:
                raise ValueError(f"Concurrent jobs disagree about output {key}")
            if task is None or task.token.cancelled:
                task = _Task(key, signature, owner, fn, max(1, int(threads)))
                self._tasks[key] = task
                self._queues.setdefault(owner, deque()).append(task)
            task.subscribers[future] = lambda: None
            if token is not None:
                detach = token.add_callback(lambda: self._cancel(task, future))
                if future in task.subscribers:
                    task.subscribers[future] = detach
                else:
                    detach()
            future.add_done_callback(
                lambda f: self._cancel(task, f) if f.cancelled() else None)
            self._pump()
        return future

    def _cancel(self, task, future):
        with self._lock:
            detach = task.subscribers.pop(future, None)
            if detach is None:
                return
            detach()
            if not future.done():
                future.set_exception(CancelledError())
            if not task.subscribers:
                task.token.cancel()
                if not task.running and self._tasks.get(task.key) is task:
                    del self._tasks[task.key]
            self._pump()

    def _pump(self):
        # Called with the lock held. Rotate between workflow queues; retain FIFO
        # within each workflow. Batches cap outstanding submissions upstream.
        if self._closed:
            return
        misses = 0
        while self._queues and len(self._running) < self._capacity:
            owner, queue = self._queues.popitem(last=False)
            while queue and not queue[0].subscribers:
                queue.popleft()
            if not queue:
                continue
            task = queue[0]
            slots = min(task.threads, self._budget)
            busy_output = any(t.key == task.key for t in self._running)
            if self._used + slots > self._budget or busy_output:
                self._queues[owner] = queue
                misses += 1
                if misses >= len(self._queues):
                    break
                continue
            misses = 0
            queue.popleft()
            if queue:
                self._queues[owner] = queue
            task.running = True
            task.allocated = slots
            self._running.add(task)
            self._used += slots
            worker = self._executor.submit(task.fn, task.token, slots)
            worker.add_done_callback(lambda f, t=task: self._finished(t, f))

    def _finished(self, task, worker):
        with self._lock:
            self._running.remove(task)
            self._used -= task.allocated
            if self._tasks.get(task.key) is task:
                del self._tasks[task.key]
            subscribers = list(task.subscribers.items())
            task.subscribers.clear()
            for future, detach in subscribers:
                detach()
                if future.done():
                    continue
                try:
                    future.set_result(worker.result())
                except BaseException as exc:
                    future.set_exception(exc)
            self._pump()

    def shutdown(self, wait=True, timeout=None):
        """Refuse new jobs, cancel every subscriber (which kills the dispatch children of
        unshared running jobs) and, when ``wait``, join the workers for at most ``timeout``
        seconds (unbounded when ``None``). Returns the output keys of jobs still running
        when it gives up: a native call that ignores its token cannot be interrupted, and
        its thread keeps running after this returns."""
        with self._lock:
            self._closed = True
            tasks = set(self._tasks.values()) | self._running
            for task in tasks:
                for future in list(task.subscribers):
                    self._cancel(task, future)
            self._queues.clear()
        self._executor.shutdown(wait=False, cancel_futures=True)
        if wait:
            join_workers(self._executor, timeout)
        with self._lock:
            return tuple(task.key for task in self._running)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.shutdown()


def join_workers(pool, timeout=None):
    """Join a shut-down ``ThreadPoolExecutor``'s worker threads for at most ``timeout``
    seconds in total; return the ones still alive. The standard library offers only an
    unbounded wait, so this reads the pool's worker set (``_threads``, stable since 3.8)."""
    deadline = None if timeout is None else time.monotonic() + timeout
    threads = list(getattr(pool, "_threads", ()))
    for thread in threads:
        thread.join(None if deadline is None else max(0.0, deadline - time.monotonic()))
    return [thread for thread in threads if thread.is_alive()]


_default = None
_default_lock = threading.Lock()
_current = ContextVar("tumult_execution", default=None)
_observer = ContextVar("tumult_run_observer", default=None)

# Process-wide coordinator pools. They are outside ExecutionService because their
# tasks may wait for leaf dispatch jobs; putting them on the leaf pool can deadlock.
IO = ThreadPoolExecutor(max_workers=4, thread_name_prefix="tumult-io")
WORKFLOWS = ThreadPoolExecutor(max_workers=4, thread_name_prefix="tumult-workflow")
HOP = ThreadPoolExecutor(max_workers=1, thread_name_prefix="tumult-hop")

# HOP keeps process-wide global state and refuses a second live instance, so every
# caller of the native library — the strength bands, and anything else that solves an
# EOS in this process — serialises on this one lock rather than on a private one.
HOP_LOCK = threading.Lock()


def get_runtime():
    """One service per server process, independent of Bokeh document imports."""
    if _current.get() is not None:
        return _current.get()
    global _default
    with _default_lock:
        if _default is None:
            _default = ExecutionService()
        return _default


@contextmanager
def execution_scope(max_threads=None):
    """Explicit lifetime for headless compute_spectrum/compute_spectra callers."""
    with ExecutionService(max_threads) as runtime:
        reset = _current.set(runtime)
        try:
            yield runtime
        finally:
            _current.reset(reset)


@contextmanager
def observe_runs(callback):
    """Report every native run folder scheduled inside this block to ``callback``.

    Batches call it from the coordinator thread before submitting the run, cache hits
    included, so a durable workflow record can name its runs while they are in flight
    instead of guessing from timestamps.
    """
    reset = _observer.set(callback)
    try:
        yield
    finally:
        _observer.reset(reset)


def note_run(output):
    observer = _observer.get()
    if observer is not None:
        observer(output)


# How long a server shutdown waits for running work before reporting what is still
# running. Cancellable work (a dispatch child, a token-checking coordinator) settles well
# within this; a native HOP call cannot be interrupted and is reported instead.
SHUTDOWN_TIMEOUT_S = 10.0


def shutdown(timeout=SHUTDOWN_TIMEOUT_S):
    """Shut down the process-wide pools within a bounded wait.

    In order: the WORKFLOWS pool stops admitting coordinators and its queued futures are
    cancelled; every live durable workflow is cancelled through its handle (queued ones
    settle as ``cancelled`` without a worker, running ones through their token); the leaf
    service cancels its jobs, which kills their dispatch children; the HOP and IO pools
    refuse new work. Then the workers of all pools are joined for at most ``timeout``
    seconds in total. The returned report says whether everything ``stopped`` and names
    the ``running`` work otherwise; a thread that ignores its token (an in-process HOP
    calculation) keeps running after this returns and is never claimed stopped.
    """
    from .submission import cancel_live_workflows
    started = time.monotonic()
    WORKFLOWS.shutdown(wait=False)               # stop admission before cancelling
    cancelled = cancel_live_workflows()
    WORKFLOWS.shutdown(wait=False, cancel_futures=True)
    for pool in (HOP, IO):
        pool.shutdown(wait=False, cancel_futures=True)
    leaf = get_runtime().shutdown(timeout=timeout)
    running = list(leaf)
    for name, pool in (("workflows", WORKFLOWS), ("hop", HOP), ("io", IO)):
        remaining = None if timeout is None else max(0.0, timeout - (time.monotonic() - started))
        running += [f"{name}:{thread.name}" for thread in join_workers(pool, remaining)]
    report = {"stopped": not running, "running": running, "cancelled_workflows": cancelled,
              "elapsed": time.monotonic() - started}
    if running:
        logging.getLogger(__name__).warning(
            "shutdown: %d thread(s) still running after %.1fs: %s",
            len(running), report["elapsed"], ", ".join(running))
    return report
