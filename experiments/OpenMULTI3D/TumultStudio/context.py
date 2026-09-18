"""AppContext - the one shared services object for the TumultStudio GUI.

``app.build`` creates the context first and composes the document's services onto it as
they come into existence: the :class:`~.execution.delivery.DeliveryService` installs the
submit entry points, the :class:`~.execution.workflows.WorkflowService` the durable-
workflow ones, and the :class:`~.workspace.persistence.WorkspacePersistence` the
workspace ones. Every tool is then built against the context rather than a bundle of
closures, and cross-tool reactions go through its event bus:

* ``ctx.bridges`` - the shared :class:`htmlbridge.BridgeSet`; each tool creates and
  owns its hidden server->client bridges here (they are all collected into the single
  document-root column by ``app.build``).
* ``ctx.tools``   - the :class:`toolmanager.ToolManager` (mutual exclusion).
* ``ctx.canvas``  - the shared :class:`canvas.SpectrumCanvas` (set right after it is
  constructed; tools are built after that).
* ``ctx.annot``   - the shared :class:`lineannot.LineAnnotations` overlay (line positions
  drawn over the spectrum); a tool attaches to it while its view is open.
* ``ctx.strengths`` - the shared :class:`tumult.linestrength.LineStrengths` provider:
  which line list, and whose strengths, the overlay ranks by and the line-list panel shows.
* ``ctx.notify``  - the footer status message (installed once the footer bridge exists;
  safe to call before that - it no-ops).
* ``ctx.run_busy`` / ``ctx.schedule`` / ``ctx.submit`` - deferred-work helpers
  (``submit`` is :meth:`~.execution.delivery.DeliveryService.submit`: a document-owned
  job on the CPU-budgeted pool, cancelled when the session closes unless admitted as
  detached).
* ``ctx.submit_io`` - the same for work that waits on the network or disk (a grid
  download), on the IO pool rather than the CPU-budgeted one.
* ``ctx.submit_operation`` - queue one public operation of :mod:`tumult.synthesis.tasks`
  on a request (a synthesis, a fit, an abundance fit, a balance, NLTE corrections) as a
  durable workflow, deliver its result to this document when it settles and stream its
  progress snapshots back; the workflow keeps running when the session closes.
  ``ctx.task_importers`` and ``ctx.task_progress`` are where a tool registers how a
  record of its operation that another document submitted is delivered to it. ``ctx.synth_modal`` is the
  :class:`synthmodal.SynthModal` that collects synthesis and fit requests.
* ``ctx.jobs``    - the :class:`jobs.JobTable`: every background operation this document
  admitted, with its queued/running/cancelling/terminal state (the Processes panel).
* ``ctx.workflow_refs`` - the :class:`jobs.WorkflowRefs` this document is attached to,
  saved beside the workspace in a slot; ``ctx.recovery`` the :class:`jobs.Recovery`
  records the recovery modal offers. ``ctx.discover_workflows``, ``ctx.attach_workflow``
  and ``ctx.import_workflow`` are the recovery actions.
* ``ctx.generation`` - the workspace generation; work captures it when admitted and a
  delivery for another generation is dropped. ``advance_generation`` marks a workspace
  boundary (a session load).
* ``ctx.restore_workspace`` - load a saved workspace slot as one document operation
  (:mod:`workspace.persistence`; ``session.py`` holds the schema and the load policy).
  ``ctx.workspace_changed`` tells that owner the live workspace moved, which is what
  keeps the autosave slot current.
* ``ctx.busy``    - raise/lower the footer progress bar without deferring anything, for
  work that already owns a tick of its own (``run_busy`` is this plus the deferral).
* ``ctx.removal`` - the :class:`removal.Removal` every sidebar remove button requests
  through: it confirms a removal that would cost something and offers Undo.
* ``ctx.execution`` - the shared :class:`tumult.synthesis.runtime.ExecutionService`
  behind the dispatch runs.
* ``ctx.emit`` / ``ctx.subscribe`` - a minimal synchronous event bus. Events in use:

  - ``"changed"``       - the workspace changed (an object was added/removed/selected,
                          units flipped, ...). Emitted from ``update_all``; tools keep
                          themselves in sync from it instead of being hand-called.
  - ``"filter-state"``  - the canvas filter existence/visibility changed
                          (``exists=bool, visible=bool``); emitted by the filter tool.
  - ``"line-picked"``   - one of the line-position indicators was clicked
                          (``lam0_aa=float, ion=str``); emitted by ``lineannot``, and what
                          opens the line-list panel at that transition's row.
"""
from __future__ import annotations

import logging

from .execution.jobs import JobTable, Recovery, WorkflowRefs


class AppContext:
    def __init__(self, doc, bridges, tools):
        self.doc = doc
        self.closed = False
        self.bridges = bridges
        self.tools = tools
        self.canvas = None
        self.annot = None      # LineAnnotations, installed by main once the canvas exists
        self.strengths = None  # LineStrengths, installed by main alongside it
        self.jobs = JobTable()
        self.workflow_refs = WorkflowRefs()
        self.recovery = Recovery()
        self.generation = object()
        self.notify = lambda *a, **k: None
        self.run_busy = lambda fn: fn()
        self.busy = lambda on: None
        self.workspace_changed = lambda: None  # dirty notification; by persistence
        self.submit = None
        self.submit_io = None
        self.submit_operation = None  # (name, operation, request, on_done, ...) -> Job; by workflows
        # a tool's delivery of a record of its kind that this document did not submit
        # (a reconnect, a loaded slot): operation name -> fn(result, handle) for the finished
        # result and fn(snapshot, handle) for an intermediate one; tools register both
        self.task_importers = {}
        self.task_progress = {}
        self.restore_workspace = None  # (snapshot, name, workflows) -> LoadReport; by persistence
        self.discover_workflows = None  # () -> offered ids; by workflows
        self.attach_workflow = None     # (id) -> Job | None; by workflows
        self.import_workflow = None     # (id) -> future | None; by workflows
        self.execution = None       # shared dispatch service; bound by main
        self.removal = None         # removal.Removal; by main
        self._subscribers = {}

    def schedule(self, fn):
        """Marshal ``fn`` onto the document thread (next tick); nothing once closed."""
        if not self.closed:
            self.doc.add_next_tick_callback(fn)

    def advance_generation(self):
        """Start a new workspace generation: outcomes of work admitted before this call
        no longer reach the workspace (they still settle their Processes rows)."""
        self.generation = object()
        return self.generation

    # --- event bus ----------------------------------------------------------
    def subscribe(self, event, fn):
        """Call ``fn(**data)`` whenever ``event`` is emitted (registration order)."""
        self._subscribers.setdefault(event, []).append(fn)
        return fn

    def emit(self, event, **data):
        """Synchronously notify the subscribers of ``event``; a failing handler is
        logged and skipped so one broken listener cannot wedge the rest of the UI."""
        for fn in self._subscribers.get(event, []):
            try:
                fn(**data)
            except Exception:
                logging.exception("handler for %r failed", event)
