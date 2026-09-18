"""Admission and delivery of one document's background work.

Every background operation a document starts is admitted here: a :class:`jobs.Job` is
created under the current workspace generation, a Processes row appears once the work
has lasted long enough to be worth showing, the work runs on one of the shared pools and
its outcome is applied on the document thread on the next tick.

Two kinds of work are distinguished. *Disposable local work* — a file load, a broadening
fit, a band of line strengths, an autosave write — runs on a pool of this server process
and is cancelled when the session closes. *Durable operations* are the dispatch
calculations: :mod:`.workflows` submits them to the workflow daemon and admits them here
as detached jobs, so closing the tab leaves them running and only their Processes row
cancels them, through the workflow handle.

Delivery is guarded twice. A closed document still settles its row (the server log and a
later reconnect see what happened) but nothing is applied to it, and an outcome admitted
under a workspace generation that has since been replaced is recorded and dropped. What
a result *means* is the caller's ``on_done``; nothing here interprets one.
"""
from __future__ import annotations

import itertools
import logging
import time

from tumult.synthesis import workflows as synth
from tumult.synthesis.runtime import IO as IO_POOL, WORKFLOWS as WORKFLOW_POOL

from ..profiling import profiled
from . import processes
from .jobs import CANCELLED, FAILED, FINISHED, Job, QUEUED, RUNNING

# Job failures are logged here from the worker thread, traceback included, whether or not
# the document that admitted them still exists; the user-facing report goes through the
# owning document's Processes row and footer alone. No root-logger handler is installed
# per document, so one session's errors never reach another session's footer.
log = logging.getLogger("tumultstudio.jobs")


class DeliveryService:
    """The document's admission and settlement of background work.

    ``on_status`` refreshes the job and loading surfaces of the page (the footer, the
    empty state, the picker spinners) without invalidating scientific tool state; it is
    the document's, not this service's, since what is on the page is the document's
    business. ``panel`` is the Processes panel this service's rows appear in.
    """

    def __init__(self, ctx, panel, *, on_status):
        self.ctx = ctx
        self.jobs = ctx.jobs
        self.panel = panel
        self._status = on_status
        self._nonce = itertools.count()
        ctx.submit = self.submit
        ctx.submit_io = self.submit_io

    # --- admission ------------------------------------------------------------------
    def admit(self, name, *, heavy=False, detached=False, path=None, delivers=True,
              operation=None):
        """Create the job row for one operation, under the current generation."""
        job = Job(next(self._nonce), name, synth.CancelToken(), self.ctx.generation,
                  heavy=heavy, detached=detached, path=str(path) if path else None,
                  delivers=delivers, operation=operation)
        self.jobs.add(job)
        # the row appears once the job has lasted SHOW_AFTER_S (see processes_html)
        try:
            self.ctx.doc.add_timeout_callback(
                lambda: None if job.terminal or self.ctx.closed else self.panel.repaint(),
                int(processes.SHOW_AFTER_S * 1000) + 20)
        except Exception:
            pass                     # the document went away
        return job

    def admission_failed(self, job, exc):
        """The pool refused the job: the same terminal path as a failure in flight, so
        the spinner clears and the row keeps the reason."""
        log.exception("%s (job %d) could not be admitted", job.name, job.pid)
        self.jobs.settle(job.pid, FAILED, error=f"{type(exc).__name__}: {exc}")
        self.ctx.notify(f"{job.name} could not start: {exc} — see Processes", level="danger")
        self._status()

    # --- disposable local work ------------------------------------------------------
    def load_async(self, name, loader, on_done, executor=IO_POOL, path=None,
                   detached=False, delivers=True):
        """Run ``loader(token)`` on ``executor`` and hand its value to ``on_done`` on the
        document thread. Returns the pool future, or ``None`` when the pool refused it.
        ``delivers`` False marks work whose outcome stays off the canvas (a grid
        download), so the empty state does not wait for it."""
        job = self.admit(name, heavy=executor is WORKFLOW_POOL, detached=detached,
                         path=path, delivers=delivers)
        token = job.token                # the loader threads it into synth.* for hard cancel
        self._status()                   # show the queued row immediately

        def work():                      # runs in a worker thread
            if token.cancelled:          # cancelled before the worker even started
                self.settle(job, None, None, None)
                return
            if job.state == QUEUED:
                job.state, job.started = RUNNING, time.time()
                self.ctx.schedule(self._status)
            try:
                result = profiled(name, lambda: loader(token))
            except synth.CancelledError:
                self.settle(job, None, None, None)  # token.cancelled -> dropped
            except Exception as exc:
                # the full traceback belongs to the server log whether or not the
                # document that asked still exists; the row gets the concise message
                log.exception("%s (job %d) failed", name, job.pid)
                self.settle(job, None, None, exc)
            else:
                self.settle(job, result, on_done, None)

        try:
            job.future = executor.submit(work)
        except Exception as exc:
            self.admission_failed(job, exc)
            return None
        return job.future

    def submit(self, name, loader, on_done, detached=False):
        """Document-owned work on the CPU-budgeted pool (a line measurement, a broadening
        fit, a normalisation): cancelled when the session closes. Dispatch calculations
        do not come through here; they are durable workflows (see :mod:`.workflows`)."""
        return self.load_async(name, loader, on_done, executor=WORKFLOW_POOL,
                               detached=detached)

    def submit_io(self, name, loader, on_done):
        """Work that waits on the network or disk (a grid download, an autosave write):
        the IO pool rather than the CPU budget, and its outcome is not a canvas object."""
        return self.load_async(name, loader, on_done, delivers=False)

    # --- settlement -----------------------------------------------------------------
    def settle(self, job, result, on_done, exc, on_error=None):
        """Deliver a job's outcome on the document thread. Called from the worker (or,
        for a queued cancellation, the cancelling thread); the state and any delivery
        are applied on the document, guarded by its lifetime and the workspace
        generation the job was admitted under. A closed document only records."""
        ctx = self.ctx

        def apply():
            if job.token.cancelled:
                self.jobs.settle(job.pid, CANCELLED)   # user cancelled: drop what arrived
            elif exc is not None:
                logs = job.handle.logs if job.handle is not None else ()
                self.jobs.settle(job.pid, FAILED, error=f"{type(exc).__name__}: {exc}",
                                 logs=logs)
            else:
                self.jobs.settle(job.pid, FINISHED)
            if ctx.closed:
                return                       # recorded; there is no document to tell
            if job.state == FAILED:
                ctx.notify(f"{job.name} failed: {exc} — see Processes", level="danger")
                if on_error is not None and job.generation is ctx.generation:
                    try:
                        on_error(exc)
                    except Exception:
                        log.exception("%s (job %d): failure not applied", job.name, job.pid)
            elif job.state == FINISHED:
                if job.generation is not ctx.generation:
                    log.info("%s (job %d) finished after its workspace was replaced; dropped",
                             job.name, job.pid)
                elif on_done is not None:
                    try:
                        on_done(result)      # add_object -> on_change -> update_all
                    except Exception as err:
                        log.exception("Failed to add %s (job %d)", job.name, job.pid)
                        ctx.notify(f"{job.name} could not be added: {err}", level="danger")
            self._status()                   # clear the spinner without broad invalidation

        if ctx.closed:
            apply()                      # no document to deliver to: record only
            return
        try:
            ctx.doc.add_next_tick_callback(apply)
        except Exception:
            apply()                      # the document went away under us

    # --- row actions ----------------------------------------------------------------
    def cancel_pending(self, pid):
        """Cancel a job from its row in the Processes panel. A still-queued job settles
        at once; a running one is told through its token (a dispatch child is killed, a
        shared dispatch job continues for its other subscribers) and its row shows
        the request until the worker settles. Its result, if any still arrives, is
        dropped in :meth:`settle`."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        job = self.jobs.request_cancel(pid)
        if job is None:
            return
        self._status()
        if job.terminal:
            self.ctx.notify(f"Cancelled {job.name}", level="info")
        else:
            self.ctx.notify(f"Cancelling {job.name}…", level="info")

    def dismiss_job(self, pid):
        """Drop a retained (failed) row from the Processes panel."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        if self.jobs.dismiss(pid) is not None:
            self._status()

    def cancel_document_work(self):
        """Closing a tab cancels this document's work (loads, fits, analysis jobs),
        except the jobs that asked not to be: every dispatch calculation is a durable
        workflow that keeps running and writes its record and result under its
        ``_workflows/`` folder, where a later document can reopen it by id."""
        for job in self.jobs.active:
            if not job.detached:
                self.jobs.request_cancel(job.pid)
