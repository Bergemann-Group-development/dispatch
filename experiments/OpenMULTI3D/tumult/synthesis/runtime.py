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
import os
import threading

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


class ExecutionService:
    """Persistent supervisors with a shared CPU budget and per-job cancellation.

    ``fn(token, threads)`` executes outside the scheduler lock. Requests with the
    same key and signature share a result; callers must treat it as read-only.
    Lowering the budget lets existing jobs drain before admitting more work.
    """

    def __init__(self, max_threads=None):
        self._capacity = max(1, os.cpu_count() or 1, int(max_threads or 1))
        self._budget = max(1, int(max_threads or self._capacity))
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

    def shutdown(self, wait=True):
        with self._lock:
            self._closed = True
            tasks = set(self._tasks.values()) | self._running
            for task in tasks:
                for future in list(task.subscribers):
                    self._cancel(task, future)
            self._queues.clear()
        self._executor.shutdown(wait=wait, cancel_futures=True)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.shutdown()


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


def shutdown():
    """Shut down the process-wide leaf service and coordinator pools."""
    get_runtime().shutdown()
    for pool in (WORKFLOWS, HOP, IO):
        pool.shutdown(wait=True, cancel_futures=True)
