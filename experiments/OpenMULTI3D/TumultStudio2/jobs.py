"""Document-owned background work: the job table behind the Processes panel, and the
result jobs that stream progress into a displayed result.

A :class:`Job` is one operation this document admitted (a file load, a synthesis, a fit,
a durable workflow) and the states it passes through::

    queued -> running -> finished | failed | cancelled
    queued -> cancelled                 (its future cancelled before a worker took it)
    queued | running -> cancelling -> cancelled

``cancelling`` is a request the worker has not yet honoured; a native call that ignores
its token stays there until it returns, and the row says so instead of claiming the work
has stopped. Failed jobs are retained with their error and native logs until dismissed;
finished and cancelled ones leave the table when they settle.
"""
from copy import deepcopy
from dataclasses import dataclass, field
import logging
from threading import Lock
import time

from tumult.synthesis.cancel import CancelToken, CancelledError

QUEUED, RUNNING, CANCELLING = "queued", "running", "cancelling"
FINISHED, FAILED, CANCELLED = "finished", "failed", "cancelled"
TERMINAL = frozenset({FINISHED, FAILED, CANCELLED})


@dataclass(eq=False)
class Job:
    """One background operation, from admission to its terminal state.

    ``generation`` is the workspace generation the job was admitted under; a delivery
    for another generation is dropped. ``future`` is the pool future when the work runs
    on a plain pool; ``handle`` the :class:`WorkflowHandle` of a durable synthesis, whose
    stable id names the on-disk record a later document can reopen."""
    pid: int
    name: str
    token: CancelToken
    generation: object
    heavy: bool = False          # confirm before cancelling (a synthesis, a fit)
    detached: bool = False       # survives the session closing; see session_closed
    path: str | None = None      # the file a load came from (picker row spinners)
    future: object = None
    handle: object = None
    state: str = QUEUED
    error: str | None = None
    logs: tuple = ()
    submitted: float = field(default_factory=time.time)
    ended: float | None = None

    @property
    def terminal(self):
        return self.state in TERMINAL

    @property
    def workflow_id(self):
        return self.handle.id if self.handle is not None else None

    def display_state(self):
        """The state to show: a queued durable workflow reads whether its worker has
        started from its record, since the pool gives no start notification."""
        if self.state == QUEUED and self.handle is not None:
            try:
                if self.handle.status == RUNNING:
                    self.state = RUNNING
            except OSError:
                pass
        return self.state


class JobTable:
    """The jobs one document admitted, in admission order (the Processes panel)."""

    def __init__(self):
        self._jobs = {}

    def add(self, job):
        self._jobs[job.pid] = job
        return job

    def get(self, pid):
        return self._jobs.get(pid)

    def __iter__(self):
        return iter(list(self._jobs.values()))

    def __len__(self):
        return len(self._jobs)

    @property
    def active(self):
        """Jobs still owed an outcome: queued, running or cancelling."""
        return [job for job in self._jobs.values() if not job.terminal]

    @property
    def retained(self):
        """Terminal jobs kept on display (failures) until dismissed."""
        return [job for job in self._jobs.values() if job.terminal]

    def paths(self):
        return [job.path for job in self.active if job.path is not None]

    def request_cancel(self, pid):
        """Ask a job to stop. A job whose future is still queued settles as cancelled
        at once; a running one becomes ``cancelling`` until its worker settles it. An
        explicit request overrides the detach opt-out. Returns the job, or ``None``
        when there is nothing to cancel."""
        job = self._jobs.get(pid)
        if job is None or job.terminal or job.state == CANCELLING:
            return None
        job.detached = False
        job.token.cancel()          # kills a dispatch child; forwards to a workflow handle
        if job.future is not None and job.future.cancel():
            self.settle(pid, CANCELLED)
        else:
            job.state = CANCELLING
        return job

    def settle(self, pid, state, error=None, logs=()):
        """Record a terminal state; failed jobs stay listed, the others leave."""
        job = self._jobs.get(pid)
        if job is None:
            return None
        job.state, job.error, job.logs, job.ended = state, error, tuple(logs), time.time()
        if state != FAILED:
            self._jobs.pop(pid, None)
        return job

    def dismiss(self, pid):
        """Drop a terminal entry from display; an active job is left alone."""
        job = self._jobs.get(pid)
        if job is not None and job.terminal:
            del self._jobs[pid]
        return job


@dataclass
class _Outcome:
    result: object = None
    error: Exception | None = None


class ResultJobs:
    """One current submission per result; workers never mutate the displayed object."""

    def __init__(self, ctx, find=None, render=None):
        # ``find``/``render`` are only consulted for result-bound jobs; a job that
        # produces a new object (``result=None``) needs neither.
        self.ctx, self.find, self.render = ctx, find, render
        self.current = {}

    def invalidate(self, uid):
        self.current.pop(uid, None)

    def submit(self, name, work, done, *, result=None, on_error=None, detached=False):
        """Submit one job. ``detached`` keeps it running when the browser tab closes,
        for work too long to babysit; its result still only reaches a live document."""
        key = result.uid if result is not None else name
        generation = object()
        self.current[key] = generation
        workspace = getattr(self.ctx, "generation", None)   # captured at admission
        working = deepcopy(result)
        lock = Lock()
        pending = []
        finished = False

        def valid():
            return (not getattr(self.ctx, "closed", False)
                    and getattr(self.ctx, "generation", None) is workspace
                    and self.current.get(key) is generation
                    and (result is None or self.find(key) is result))

        def apply(snapshot):
            result.__dict__.update(snapshot.__dict__)
            self.render(result)

        def flush():
            with lock:
                snapshot = pending.pop() if pending else None
            if snapshot is not None and valid() and not finished:
                apply(snapshot)

        def progress(value):
            # Coalesce before copying; published arrays/spectra are independent of the worker.
            with lock:
                if pending:
                    return
                pending.append(deepcopy(value))
            self.ctx.schedule(flush)

        def stopped(status):
            nonlocal finished
            if finished or not valid():
                return
            finished = True
            self.current.pop(key, None)
            if result is not None:
                result.termination = status
                if hasattr(result, "converged"):
                    result.converged = False
                self.render(result)

        def run(token):
            token = token or CancelToken()
            detach = token.add_callback(lambda: self.ctx.schedule(lambda: stopped("cancelled")))
            try:
                token.check()
                return _Outcome(work(working, token, progress))
            except CancelledError:
                raise
            except Exception as exc:
                # the toast shows only the message; the traceback goes to the server log
                logging.exception("%s failed", name)
                return _Outcome(error=exc)
            finally:
                detach()

        def complete(outcome):
            nonlocal finished
            if not valid() or finished:
                return
            if outcome.error is not None:
                stopped("failed")
                if on_error is not None:
                    on_error(outcome.error)      # the submitter has its own error surface
                else:
                    self.ctx.notify(f"{name} failed: {outcome.error}")
                return
            finished = True
            self.current.pop(key, None)
            if result is not None:
                apply(outcome.result)
                done(result)
            else:
                done(outcome.result)

        future = self.ctx.submit(name, run, complete, detached=detached)
        if future is not None:
            future.add_done_callback(
                lambda f: self.ctx.schedule(lambda: stopped("cancelled")) if f.cancelled() else None)
        return future
