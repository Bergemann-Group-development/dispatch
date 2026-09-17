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

A durable workflow (every dispatch calculation: a synthesis, a fit, a tool run) runs in
the workflow daemon and outlives the document: :class:`WorkflowRefs` remembers which ones
this document is attached to (saved beside the workspace in a slot), and :class:`Recovery`
lists the settled records the recovery modal offers to import or forget. The daemon's
running and queued workflows that this document does not follow, and dispatch runs
started outside any workflow (:class:`ExternalRun`), are shown beside its own jobs. Nothing
here resumes a native calculation; a record whose daemon exited is shown as interrupted.
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
    on a plain pool; ``handle`` the :class:`WorkflowHandle` of a durable workflow, whose
    stable id names the on-disk record a later document can reopen. ``threads`` is what
    the job runs on: one for a Python job of this server, the daemon's allocation for a
    workflow. ``started`` is when it began running (epoch seconds)."""
    pid: int
    name: str
    token: CancelToken
    generation: object
    heavy: bool = False          # confirm before cancelling (a synthesis, a fit)
    detached: bool = False       # survives the session closing; see session_closed
    path: str | None = None      # the file a load came from (picker row spinners)
    delivers: bool = True        # its outcome lands on the canvas (a load, a synthesis)
    operation: str | None = None  # the task kind of a durable workflow
    future: object = None
    handle: object = None
    state: str = QUEUED
    error: str | None = None
    logs: tuple = ()
    submitted: float = field(default_factory=time.time)
    started: float | None = None
    ended: float | None = None
    threads: int = 1

    @property
    def terminal(self):
        return self.state in TERMINAL

    @property
    def workflow_id(self):
        return self.handle.id if self.handle is not None else None

    def display_state(self, summary=None):
        """The state to show. A queued durable workflow is running once the daemon says
        so: its ``summary`` from the daemon's live jobs when given, else its record."""
        if self.handle is not None and summary is not None:
            self.threads = summary.get("threads") or self.threads
            self.started = summary.get("started") or self.started
        if self.state == QUEUED and self.handle is not None:
            try:
                running = (summary.get("state") == RUNNING if summary is not None
                           else self.handle.status == RUNNING)
            except OSError:
                running = False
            if running:
                self.state = RUNNING
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

    def find_workflow(self, workflow_id):
        """The active job attached to workflow ``workflow_id``, if any."""
        return next((job for job in self.active if job.workflow_id == workflow_id), None)

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
class WorkflowRef:
    """A durable workflow this document is attached to: enough to reopen it later
    (``folder``), name it, and know whether its result already reached the workspace."""
    id: str
    folder: str
    name: str
    imported: bool = False

    def to_dict(self):
        return {"id": self.id, "folder": self.folder, "name": self.name,
                "imported": self.imported}

    @classmethod
    def from_dict(cls, value):
        return cls(str(value["id"]), str(value["folder"]), str(value["name"]),
                   bool(value.get("imported", False)))


class WorkflowRefs:
    """The workflow references kept beside the workspace: every workflow this document
    submitted or attached to, and which of them have delivered. A save slot stores them
    (``store.save_session(..., workflows=refs.to_list())``); a load replaces them."""

    def __init__(self):
        self._refs = {}

    def attach(self, handle, name, imported=False):
        ref = self._refs.get(handle.id)
        if ref is None:
            ref = self._refs[handle.id] = WorkflowRef(handle.id, str(handle.folder), name)
        ref.imported = ref.imported or imported
        return ref

    def mark_imported(self, workflow_id):
        ref = self._refs.get(workflow_id)
        if ref is not None:
            ref.imported = True
        return ref

    def imported(self, workflow_id):
        ref = self._refs.get(workflow_id)
        return ref is not None and ref.imported

    def get(self, workflow_id):
        return self._refs.get(workflow_id)

    def __contains__(self, workflow_id):
        return workflow_id in self._refs

    def __iter__(self):
        return iter(list(self._refs.values()))

    def __len__(self):
        return len(self._refs)

    def replace(self, refs):
        """Adopt a slot's references (dicts or :class:`WorkflowRef`) wholesale."""
        self._refs = {}
        for ref in refs:
            ref = ref if isinstance(ref, WorkflowRef) else WorkflowRef.from_dict(ref)
            self._refs[ref.id] = ref

    def to_list(self):
        return [ref.to_dict() for ref in self._refs.values()]


@dataclass(eq=False)
class Recoverable:
    """An on-disk workflow record offered for recovery, as read once when listed: a live
    one can be attached (a new generation-bound subscription with cancellation), a
    finished one imported, and any of them forgotten. ``interrupted`` is a record left
    queued or running by a process that no longer exists; the only honest offer is its
    logs. ``referenced`` says a loaded slot named it."""
    handle: object
    name: str
    state: str
    referenced: bool = False
    error: str | None = None
    logs: tuple = ()
    phase: str | None = None
    elapsed: float = 0.0
    runs: int = 0
    kind: str | None = None

    @property
    def id(self):
        return self.handle.id

    @property
    def live(self):
        return self.state in (QUEUED, RUNNING)

    @classmethod
    def from_handle(cls, handle, referenced=False):
        record = handle.to_dict()
        state = handle.status
        ended = record.get("ended") or time.time()
        try:
            kind = handle.kind.name
        except KeyError:
            kind = None
        return cls(handle, handle.label, state, referenced, record.get("error"),
                   handle.logs, record.get("phase"),
                   max(0.0, ended - record.get("submitted", ended)),
                   len(record.get("runs") or ()), kind)


@dataclass(eq=False)
class ExternalRun:
    """A dispatch run found by its own record (``job.nml``, see
    ``tumult.synthesis.jobrecord``) rather than by a workflow: one started from a shell,
    a notebook or a batch script, as the workflow daemon lists it. Its row can end it,
    and a departure solve can be handed to the synthesis modal as saved departures."""
    folder: str
    name: str
    state: str          # running | submitted
    pid: int | None = None
    host: str | None = None
    slurm_id: str | None = None
    input: str | None = None
    started: float | None = None
    logs: tuple = ()
    departures: bool = False
    threads: int | None = None

    @property
    def id(self):
        from tumult.synthesis.jobrecord import run_id
        return run_id(self.folder)

    @property
    def live(self):
        return self.state in ("running", "submitted")

    @classmethod
    def from_summary(cls, summary):
        """From an ``external`` job of the daemon's list."""
        from pathlib import Path
        folder = Path(summary["folder"])
        log = folder.with_suffix(".log")
        return cls(str(folder), folder.name, summary["state"], summary.get("pid"),
                   summary.get("host"), summary.get("slurm_id"), summary.get("input"),
                   summary.get("started"), (str(log),) if log.is_file() else (),
                   (folder / "save").is_dir() and any((folder / "save").glob("pop_*")),
                   summary.get("threads"))


class Recovery:
    """The workflow records the recovery modal offers (:class:`Recoverable`), keyed by
    workflow id, in the order offered."""

    def __init__(self):
        self._entries = {}

    def offer(self, entry):
        self._entries[entry.id] = entry
        return entry

    def get(self, workflow_id):
        return self._entries.get(workflow_id)

    def forget(self, workflow_id):
        return self._entries.pop(workflow_id, None)

    def __contains__(self, workflow_id):
        return workflow_id in self._entries

    def __iter__(self):
        return iter(list(self._entries.values()))

    def __len__(self):
        return len(self._entries)


@dataclass
class _Outcome:
    result: object = None
    error: Exception | None = None


def update_kept_name(result, snapshot):
    """Update ``result`` in place from a worker's ``snapshot`` of it, except its name: the
    document names its objects (:mod:`naming`), a calculation does not rename them."""
    result.__dict__.update({k: v for k, v in snapshot.__dict__.items() if k != "name"})


class ResultJobs:
    """One current submission per result; workers never mutate the displayed object.
    A result owns at most one live durable workflow (``live``): superseding or
    invalidating the submission cancels it through the Processes path
    (:meth:`JobTable.request_cancel`), so its record settles as cancelled."""

    def __init__(self, ctx, find=None, render=None):
        # ``find``/``render`` are only consulted for result-bound jobs; a job that
        # produces a new object (``result=None``) needs neither.
        self.ctx, self.find, self.render = ctx, find, render
        self.current = {}
        self.live = {}

    def _cancel_live(self, key):
        job = self.live.pop(key, None)
        if job is not None:
            self.ctx.jobs.request_cancel(job.pid)

    def invalidate(self, uid):
        self.current.pop(uid, None)
        self._cancel_live(uid)

    def submit(self, name, work, done, *, result=None, on_error=None, detached=False):
        """Submit one job. ``detached`` keeps it running when the browser tab closes,
        for work too long to babysit; its result still only reaches a live document."""
        key = result.uid if result is not None else name
        generation = object()
        self.current[key] = generation
        self._cancel_live(key)
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
            update_kept_name(result, snapshot)
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

    def submit_operation(self, name, operation, request, done, *, kwargs=None, result=None,
                         on_progress=None, on_error=None, meta=None):
        """Submit one public ``operation`` (a name in :data:`tumult.synthesis.tasks.KINDS`)
        on ``request`` with its keyword arguments as a durable workflow, through
        ``ctx.submit_operation``. It outlives the tab; the record's progress snapshots
        reach ``result`` in place (or ``on_progress`` when the job creates its result), the
        finished result reaches ``done`` and a failure ``on_error`` (the Processes row
        already says why). Returns the job, or ``None`` when submission was refused.
        Submitting again for the same ``key`` (``result.uid`` or ``name``) while an
        earlier one is still running is a deliberate rerun (changed settings, say): it
        replaces the earlier submission's ``current`` entry first, so nothing of that one
        reaches ``result`` any more, then cancels its workflow as its Processes row would.
        Its record settles as cancelled and is not offered as a result to import. A
        double *click*, rather than a deliberate rerun, is guarded upstream of here — the
        submit button disables itself and the dialog closes on the first click, before a
        second can be sent."""
        key = result.uid if result is not None else name
        generation = object()
        self.current[key] = generation       # the earlier submission is stale from here
        self._cancel_live(key)
        workspace = getattr(self.ctx, "generation", None)   # captured at admission
        finished = False
        job = None

        def release():
            if job is not None and self.live.get(key) is job:
                del self.live[key]

        def valid():
            return (not getattr(self.ctx, "closed", False)
                    and getattr(self.ctx, "generation", None) is workspace
                    and self.current.get(key) is generation
                    and (result is None or self.find(key) is result))

        def apply(snapshot):
            update_kept_name(result, snapshot)
            self.render(result)

        def progress(snapshot):          # document thread, one decoded snapshot
            if finished or not valid():
                return
            if result is not None:
                apply(snapshot)
            elif on_progress is not None:
                on_progress(snapshot)

        def stopped(status):
            nonlocal finished
            if finished or not valid():
                return
            finished = True
            self.current.pop(key, None)
            release()
            if result is not None:
                result.termination = status
                if hasattr(result, "converged"):
                    result.converged = False
                self.render(result)

        def complete(outcome):
            nonlocal finished
            if not valid() or finished:
                return
            finished = True
            self.current.pop(key, None)
            release()
            if result is not None:
                apply(outcome)
                done(result)
            else:
                done(outcome)

        def failed(exc):
            if finished or not valid():
                return
            stopped("failed")
            if on_error is not None:
                on_error(exc)

        job = self.ctx.submit_operation(name, operation, request, complete, kwargs=kwargs,
                                        meta=meta, on_progress=progress, on_error=failed)
        if job is None:
            stopped("failed")            # refused at admission; the row keeps the reason
            return None
        if not finished:
            self.live[key] = job
        job.token.add_callback(lambda: self.ctx.schedule(lambda: stopped("cancelled")))
        return job
