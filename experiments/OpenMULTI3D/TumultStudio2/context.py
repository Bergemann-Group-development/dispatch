"""AppContext - the one shared services object for the TumultStudio2 GUI.

``app.build`` creates the context first and installs the document-level services on it
as they come into existence; every tool is then built against the context rather than a
bundle of closures, and cross-tool reactions go through its event bus:

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
  (``submit`` is the background loader ``app.load_async``: a document-owned job on the
  IO pool, cancelled when the session closes unless admitted as detached).
* ``ctx.submit_io`` - the same for work that waits on the network or disk (a grid
  download), on the IO pool rather than the CPU-budgeted one.
* ``ctx.submit_workflow`` - queue one durable synthesis request (every synthesis of the
  modal) and deliver its result to this document when it settles; the workflow keeps
  running when the session closes. ``ctx.submit_task`` does the same for the tool
  calculations of :mod:`tumult.synthesis.tasks` (a fit, an abundance run, a balance,
  NLTE corrections), streaming their progress snapshots back; ``ctx.task_importers``
  and ``ctx.task_progress`` are where a tool registers how a record of its kind that
  another document submitted is delivered to it. ``ctx.synth_modal`` is the
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
  (see ``session.py`` for the schema and the load policy).
* ``ctx.busy``    - raise/lower the footer progress bar without deferring anything, for
  work that already owns a tick of its own (``run_busy`` is this plus the deferral).
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

from jobs import JobTable, Recovery, WorkflowRefs


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
        self.submit = None
        self.submit_io = None
        self.submit_workflow = None  # (name, request, on_done) -> WorkflowHandle; by main
        self.submit_task = None      # (name, kind, args, on_done, ...) -> Job; by main
        # a tool's delivery of a record of its kind that this document did not submit
        # (a reconnect, a loaded slot): kind name -> fn(result, handle) for the finished
        # result and fn(snapshot, handle) for an intermediate one; tools register both
        self.task_importers = {}
        self.task_progress = {}
        self.restore_workspace = None  # (snapshot, name, workflows) -> LoadReport; by main
        self.discover_workflows = None  # () -> offered ids; by main
        self.attach_workflow = None     # (id) -> Job | None; by main
        self.import_workflow = None     # (id) -> future | None; by main
        self.execution = None       # shared dispatch service; bound by main
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
