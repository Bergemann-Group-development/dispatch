"""One document's relationship to the durable workflows of :mod:`tumult.synthesis`.

A dispatch calculation is not this document's work. It is submitted as a durable
workflow, runs in the workflow daemon, writes its record and its result under a
``_workflows/`` folder, and outlives the tab that asked for it. What a document owns is
its *attachment* to such a workflow, which is what this module is:

* submission of one public operation as a workflow, with the Processes row, the
  cancellation path and the reference kept beside the workspace;
* the subscription that delivers a settled handle's outcome exactly once, through
  :class:`~.delivery.DeliveryService`, and marks the record imported so a later load
  does not offer the same result again;
* progress following while a workflow runs, decoded off the document thread and stopped
  by the job, the document or a workspace generation change;
* recovery after a reconnect or a workspace load: discovering settled records under the
  configured roots, reattaching to live ones, importing finished ones, and reporting
  those whose record is gone.

Nothing here schedules, cancels or stores a calculation itself; ``tumult.synthesis``
does all of that. The daemon's own jobs are read through the Processes panel's
:class:`~.processes.LiveJobs`.
"""
from __future__ import annotations

import logging
import subprocess
import sys

from tumult.synthesis import daemon as workflow_daemon
from tumult.synthesis import submission  # (open workflows on reconnect)
from tumult.synthesis import submit_operation as submit_operation_record
from tumult.synthesis import workflows as synth
from tumult.synthesis.runtime import IO as IO_POOL

from ..ui.htmlbridge import pulse
from ..workspace import session
from . import processes
from .jobs import CANCELLED, ExternalRun, FINISHED, QUEUED, RUNNING, Recoverable

log = logging.getLogger("tumultstudio.jobs")

PROGRESS_POLL_MS = 2000        # Slurm records only; the daemon reports its snapshots


class WorkflowService:
    """This document's durable-workflow attachments, recovery modal and job details.

    ``add_spectrum`` is how a synthesis result reaches the canvas — the one thing about a
    workflow's outcome this service has to be told, since every other operation's result
    goes to the importer its own tool registered on the context. ``workflow_roots``
    names where records are discovered.
    """

    def __init__(self, ctx, delivery, panel, *, add_spectrum, workflow_roots, on_status):
        self.ctx = ctx
        self.delivery = delivery
        self.panel = panel
        self.jobs = ctx.jobs
        self._add_spectrum = add_spectrum
        self._roots = workflow_roots
        self._status = on_status
        # the job details and recovery modals, filled before they open
        self._details = ctx.bridges.html("job-details", name="job-details")
        self._job_modal = ctx.bridges.dialog("job-modal")
        self._recovery = ctx.bridges.html("recovery-list", name="recovery-list")
        self._recovery_modal = ctx.bridges.dialog("recovery-modal")
        ctx.submit_operation = self.submit_operation
        ctx.discover_workflows = self.discover
        ctx.attach_workflow = self.attach
        ctx.import_workflow = self.import_result

    # --- delivery of a workflow's outcome -------------------------------------------
    def _importer(self, handle):
        """The delivery of a workflow this document did not submit: a synthesis adds
        every spectrum of its result under the name and annotations recorded at
        submission; any other operation goes to the importer its tool registered."""
        kind = handle.kind.name
        if kind != "synthesis.synthesize":
            importer = self.ctx.task_importers.get(kind)
            if importer is None:
                raise ValueError(f"no tool imports a {kind} result")
            return lambda result: importer(result, handle)

        def deliver(result):
            with self.ctx.canvas.batch_update():
                for index, spec in enumerate(result.spectra, start=1):
                    suffix = f"·{index}" if len(result.spectra) > 1 else ""
                    self._add_spectrum(spec, handle.label + suffix, meta=handle.meta)
        return deliver

    def _progress_importer(self, handle):
        """The progress delivery of a tool record this document did not submit, or
        ``None`` for a kind without progress or without a registered tool."""
        kind = handle.kind
        follow = self.ctx.task_progress.get(kind.name) if kind.progress else None
        return None if follow is None else (lambda snapshot: follow(snapshot, handle))

    def _subscribe(self, job, handle, on_done, on_error=None):
        """Deliver ``handle``'s outcome to ``job`` once it settles; the delivery marks
        the workflow imported so a later load does not offer its result again."""
        def deliver(result):
            on_done(result)
            self.ctx.workflow_refs.mark_imported(handle.id)

        def settled(handle):             # worker thread (or the canceller's), never waits
            try:
                result = handle.result()
            except synth.CancelledError:
                self.delivery.settle(job, None, None, None)
            except Exception as exc:
                self.delivery.settle(job, None, None, exc, on_error)
            else:
                self.delivery.settle(job, result, deliver, None)

        handle.add_done_callback(settled)

    def _follow_progress(self, job, handle, on_progress):
        """While ``job`` runs, hand every new progress snapshot of its record to
        ``on_progress`` on the document thread. The daemon reports each snapshot's stamp
        (a Slurm record is polled instead); the snapshot is decoded off the document
        thread, since it can carry a spectrum, one read at a time, and a stamp that
        arrives during a read is read next. Ends with the job, the document or the
        workspace generation."""
        if on_progress is None or not handle.kind.progress:
            return
        ctx = self.ctx
        state = {"stamp": None, "busy": False, "stop": None}

        def stop():
            stopper, state["stop"] = state["stop"], None
            if stopper is not None:
                try:
                    stopper()
                except Exception:
                    pass                 # already removed, or the document is gone

        def tick():
            if job.terminal or ctx.closed or job.generation is not ctx.generation:
                stop()
                return
            if state["busy"]:
                return
            try:
                stamp = handle.progress_stamp()
            except OSError:
                return
            if stamp is None or stamp == state["stamp"]:
                return
            state["busy"] = True

            def read():                  # IO pool
                try:
                    snapshot = handle.progress_result()
                except Exception:
                    log.exception("%s (job %d): progress snapshot unreadable", job.name, job.pid)
                    snapshot = None

                def deliver():
                    state["busy"], state["stamp"] = False, stamp
                    if snapshot is None or job.terminal or job.generation is not ctx.generation:
                        return
                    try:
                        on_progress(snapshot)
                    except Exception:
                        log.exception("%s (job %d): progress not applied", job.name, job.pid)
                    tick()               # a newer snapshot written during the read
                ctx.schedule(deliver)
            IO_POOL.submit(read)

        if handle.detached == "slurm":
            callback = ctx.doc.add_periodic_callback(tick, PROGRESS_POLL_MS)
            state["stop"] = lambda: ctx.doc.remove_periodic_callback(callback)
        else:
            state["stop"] = handle.on_progress(lambda stamp: ctx.schedule(tick))
            ctx.schedule(tick)           # a snapshot written before this document followed

    # --- submission ------------------------------------------------------------------
    def submit_operation(self, name, operation, request, on_done, *, kwargs=None, meta=None,
                         on_progress=None, on_error=None):
        """Queue one public operation (:mod:`tumult.synthesis.tasks`) on ``request`` as a
        durable workflow and deliver its result to ``on_done`` once its handle settles, a
        failure to ``on_error`` after the row has recorded it, and every progress snapshot
        of the record to ``on_progress`` while it runs. Submission happens here, once;
        nothing waits on the workflow from another pool thread, and the document thread
        never blocks on the result. The workflow is detached: closing the tab leaves it
        running, and only its Processes row cancels it, through the handle. ``name`` and
        ``meta`` go into the record, so a later document lists and imports it as it was
        submitted, and the workflow is referenced beside the workspace from now on. The
        driver and the record root follow the request's configuration. Returns the job,
        or ``None`` when submission failed (the row keeps the reason)."""
        job = self.delivery.admit(name, heavy=True, detached=True, operation=operation)
        try:
            handle = submit_operation_record(operation, request, label=name, meta=meta,
                                             **(kwargs or {}))
        except Exception as exc:
            self.delivery.admission_failed(job, exc)
            return None
        job.handle = handle
        # the row's cancel token forwards to the handle, so a settle can tell a cancelled
        # workflow's late outcome from a wanted one the same way it does for loads
        job.token.add_callback(handle.cancel)
        self.ctx.workflow_refs.attach(handle, name)
        self._status()
        self.ctx.workspace_changed()
        self._subscribe(job, handle, on_done, on_error)
        self._follow_progress(job, handle, on_progress)
        return job

    # --- recovery of durable workflows after a reconnect ------------------------------
    def _offer(self, handle, referenced=False):
        """List ``handle`` for recovery unless this document already follows it or has
        imported its result. Returns the entry, or ``None``."""
        if (self.jobs.find_workflow(handle.id) is not None
                or self.ctx.workflow_refs.imported(handle.id)):
            return None
        try:
            entry = Recoverable.from_handle(handle, referenced)
        except (OSError, ValueError, KeyError) as exc:
            log.warning("Workflow record %s cannot be read: %s", handle.folder, exc)
            return None
        return self.ctx.recovery.offer(entry)

    def discover(self, limit=20):
        """Ask the daemon for the newest settled records under the configured roots and
        offer the ones this document has not imported in the recovery modal, which opens
        on them. Running and queued workflows are rows of the panel already. Returns the
        ids offered."""
        try:
            workflow_daemon.client()
            view = self.panel.live.view()
            reply = view.refresh(self._roots(), limit) if view is not None else {"jobs": []}
        except (OSError, workflow_daemon.DaemonError) as exc:
            self.ctx.notify(f"The workflow daemon cannot be reached: {exc}", level="danger")
            return []
        offered = []
        for summary in reply["jobs"]:
            if summary["type"] != "workflow" or summary["state"] in (QUEUED, RUNNING):
                continue
            handle = submission.open_workflow(summary["folder"])
            entry = self._offer(handle, referenced=handle.id in self.ctx.workflow_refs)
            if entry is not None:
                offered.append(entry.id)
        self.show_recovery()
        return offered

    def show_recovery(self):
        self._render_recovery()
        self.panel.repaint()
        pulse(self._recovery_modal, "open")

    def _render_recovery(self):
        self._recovery.value = processes.recovery_html(self.ctx.recovery)

    def _reopen(self, workflow_id):
        """The handle for an offered, referenced or live workflow; ``None`` if it is
        unknown."""
        entry = self.ctx.recovery.get(workflow_id)
        if entry is not None:
            return entry.handle
        ref = self.ctx.workflow_refs.get(workflow_id)
        if ref is not None:
            return submission.open_workflow(ref.folder)
        summary = self.panel.live.get(workflow_id)
        return submission.open_workflow(summary["folder"]) if summary is not None else None

    def attach(self, workflow_id):
        """Follow a live workflow from here: a new job under the current generation whose
        row can cancel it through the handle and whose completion is delivered once. Live
        means queued or running in the daemon, or a Slurm job that still exists; a record
        without one cannot be attached and nothing is resumed."""
        handle = self._reopen(workflow_id)
        if handle is None:
            return None
        job = self.jobs.find_workflow(workflow_id)
        if job is not None:
            return job
        try:
            state = handle.status
        except OSError as exc:
            self.ctx.notify(f"Workflow {workflow_id} cannot be read: {exc}", level="danger")
            return None
        if handle.done or state not in (QUEUED, RUNNING):
            self.ctx.notify(f"{handle.label} is {state}; it cannot be attached", level="warn")
            return None
        try:
            deliver = self._importer(handle)
        except ValueError as exc:
            self.ctx.notify(f"{handle.label} cannot be attached: {exc}", level="warn")
            return None
        job = self.delivery.admit(handle.label, heavy=True, detached=True)
        job.handle = handle
        job.started = handle.to_dict().get("started")
        job.token.add_callback(handle.cancel)
        self.ctx.workflow_refs.attach(handle, handle.label)
        self.ctx.recovery.forget(workflow_id)
        pulse(self._job_modal, "close")
        self._subscribe(job, handle, deliver)
        self._follow_progress(job, handle, self._progress_importer(handle))
        self._status()
        self.ctx.workspace_changed()
        return job

    def import_result(self, workflow_id):
        """Load a finished workflow's saved result off the document thread and add its
        spectra once; a workflow already imported is refused."""
        handle = self._reopen(workflow_id)
        if handle is None:
            return None
        if self.ctx.workflow_refs.imported(workflow_id):
            self.ctx.notify(f"{handle.label} was already imported", level="warn")
            self.ctx.recovery.forget(workflow_id)
            self._status()
            return None
        try:
            deliver = self._importer(handle)
        except ValueError as exc:
            self.ctx.notify(f"{handle.label} cannot be imported: {exc}", level="warn")
            return None

        def done(result):
            deliver(result)
            self.ctx.workflow_refs.attach(handle, handle.label, imported=True)

        self.ctx.recovery.forget(workflow_id)
        self._render_recovery()
        return self.delivery.load_async(f"import {handle.label}",
                                        lambda token=None: handle.result(), done)

    def forget(self, workflow_id):
        if self.ctx.recovery.forget(workflow_id) is not None:
            self._render_recovery()
            self._status()

    def reattach(self, refs):
        """Adopt a loaded slot's workflow references: a live one is followed again (a job
        already following it is re-subscribed under the new generation), a finished one
        not yet imported is offered for import, an interrupted, failed or cancelled one is
        shown as such. Returns the outcome lines for the load report and the references
        whose record is gone."""
        self.ctx.workflow_refs.replace(refs)
        outcomes, gone = [], []
        for ref in self.ctx.workflow_refs:
            handle = submission.open_workflow(ref.folder)
            try:
                state = handle.status
            except (OSError, ValueError) as exc:
                gone.append(session.Skipped("workflow", ref.id, ref.name, ref.folder,
                                            "no record" if isinstance(exc, FileNotFoundError)
                                            else str(exc)))
                continue
            job = self.jobs.find_workflow(ref.id)
            if job is not None and not handle.done:
                job.generation = self.ctx.generation      # the slot names it: subscribe anew
                self.ctx.recovery.forget(ref.id)
                outcomes.append(f"{ref.name} ({state}, followed)")
            elif not handle.done and state in (QUEUED, RUNNING):
                self.attach(ref.id)
                outcomes.append(f"{ref.name} ({state}, reattached)")
            elif ref.imported:
                continue                              # its spectra are in the slot itself
            elif self._offer(handle, referenced=True) is not None:
                outcomes.append(f"{ref.name} ({state}"
                                + (", result to import)" if state == FINISHED else ")"))
        return outcomes, gone

    # --- rows that are not this document's jobs --------------------------------------
    def external_run(self, run_id):
        """The dispatch run ``run_id`` as the daemon lists it, or ``None``."""
        summary = self.panel.live.get(run_id)
        if summary is None or summary["type"] != "external":
            return None
        return ExternalRun.from_summary(summary)

    def cancel_run(self, run_id):
        """End an external run through the daemon, which signals it only when the process
        table confirms the recorded pid still runs that input (or scancels it)."""
        entry = self.external_run(run_id)
        if entry is None:
            return
        try:
            state = workflow_daemon.client().request("cancel", job=run_id)["state"]
        except (OSError, workflow_daemon.DaemonError) as exc:
            self.ctx.notify(f"{entry.name} could not be ended: {exc}", level="danger")
            return
        if state == CANCELLED:
            self.ctx.notify(f"Ended {entry.name}", level="info")
            pulse(self._job_modal, "close")
        else:
            self.ctx.notify(f"{entry.name} could not be signalled from here; it is {state}",
                            level="warn")
        self.panel.repaint()

    def cancel_workflow(self, workflow_id):
        """Cancel a daemon workflow this document does not follow (another tab's, a
        script's) from its row."""
        summary = self.panel.live.get(workflow_id)
        folder = summary["folder"] if summary else None
        if folder is None:
            entry = self.ctx.recovery.get(workflow_id)
            folder = entry.handle.folder if entry is not None else None
        if folder is None:
            return
        handle = submission.open_workflow(folder)
        try:
            told = handle.cancel()
        except (OSError, workflow_daemon.DaemonError) as exc:
            self.ctx.notify(f"{handle.label} could not be cancelled: {exc}", level="danger")
            return
        self.ctx.notify(f"Cancelling {handle.label}…" if told
                        else f"{handle.label} is {handle.status}",
                        level="info" if told else "warn")
        pulse(self._job_modal, "close")
        self.panel.repaint()

    # --- the job details modal --------------------------------------------------------
    def _details_target(self, uid):
        """What a details button names: ``("job", Job)``, ``("workflow", handle,
        summary, recoverable)`` or ``("run", ExternalRun)``; ``None`` if it is gone."""
        if uid.startswith("job-"):
            try:
                job = self.jobs.get(int(uid[4:]))
            except ValueError:
                return None
            return None if job is None else ("job", job)
        summary = self.panel.live.get(uid)
        if summary is not None and summary["type"] == "external":
            return ("run", ExternalRun.from_summary(summary))
        recoverable = self.ctx.recovery.get(uid)
        if summary is not None:
            return ("workflow", submission.open_workflow(summary["folder"]), summary, None)
        if recoverable is not None:
            return ("workflow", recoverable.handle, None, recoverable)
        ref = self.ctx.workflow_refs.get(uid)
        if ref is not None:
            return ("workflow", submission.open_workflow(ref.folder), None, None)
        return None

    def show_details(self, uid):
        """Fill and open the job details modal for a row's details button."""
        target = self._details_target(uid)
        if target is None:
            self.ctx.notify("That job is no longer listed", level="warn")
            self.panel.repaint()
            return
        try:
            if target[0] == "run":
                body = processes.run_details_html(target[1])
            elif target[0] == "job" and target[1].handle is None:
                body = processes.job_details_html(target[1])
            elif target[0] == "job":
                job = target[1]
                body = processes.workflow_details_html(
                    job.handle, self.panel.live.get(job.workflow_id), job=job)
            else:
                _, handle, summary, recoverable = target
                follower = self.jobs.find_workflow(handle.id)
                body = processes.workflow_details_html(handle, summary, job=follower,
                                                       recoverable=recoverable)
        except (OSError, ValueError, KeyError) as exc:
            self.ctx.notify(f"The details cannot be read: {exc}", level="danger")
            return
        self._details.value = body
        pulse(self._job_modal, "open")

    def open_folder(self, uid):
        """Show a job's folder in the file manager of the machine the server runs on."""
        target = self._details_target(uid)
        if target is None:
            return
        if target[0] == "run":
            folder = target[1].folder
        elif target[0] == "job":
            folder = target[1].handle.folder if target[1].handle is not None else None
        else:
            folder = target[1].folder
        if folder is None:
            return
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        try:
            subprocess.Popen([opener, str(folder)], stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
        except OSError as exc:
            self.ctx.notify(f"Could not open {folder}: {exc}", level="warn")
