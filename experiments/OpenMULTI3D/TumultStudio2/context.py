"""AppContext - the one shared services object for the TumultStudio2 GUI.

Before the Phase 4 consolidation, ``main.build()`` wired every tool with a bundle of
closure lambdas (AbundanceTool took 15 of them), and cross-tool reactions were
hand-threaded through ``update_all``. The context replaces that plumbing:

* ``ctx.bridges`` - the shared :class:`htmlbridge.BridgeSet`; each tool creates and
  owns its hidden server->client bridges here (they are all collected into the single
  document-root column by main).
* ``ctx.tools``   - the :class:`toolmanager.ToolManager` (mutual exclusion).
* ``ctx.canvas``  - the shared :class:`canvas.SpectrumCanvas` (set by main right after
  it is constructed; tools are built after that).
* ``ctx.annot``   - the shared :class:`lineannot.LineAnnotations` overlay (line positions
  drawn over the spectrum); a tool attaches to it while its view is open.
* ``ctx.strengths`` - the shared :class:`linestrength.LineStrengths` provider: which line
  list, and whose strengths, the overlay ranks by and the line-list panel shows.
* ``ctx.notify``  - the footer status message (installed by main once the footer
  bridge exists; safe to call before that - it no-ops).
* ``ctx.run_busy`` / ``ctx.schedule`` / ``ctx.submit`` - deferred-work helpers
  (installed by main; ``submit`` is the background loader ``load_async``).
* ``ctx.submit_workflow`` - queue one durable synthesis request and deliver its result
  to this document when it settles.
* ``ctx.jobs``    - the :class:`jobs.JobTable`: every background operation this document
  admitted, with its queued/running/cancelling/terminal state (the Processes panel).
* ``ctx.workflow_refs`` - the :class:`jobs.WorkflowRefs` this document is attached to,
  saved beside the workspace in a slot; ``ctx.recovery`` the :class:`jobs.Recovery`
  records offered after a reconnect. ``ctx.discover_workflows``, ``ctx.attach_workflow``
  and ``ctx.import_workflow`` (installed by main) are the recovery actions.
* ``ctx.generation`` - the workspace generation; work captures it when admitted and a
  delivery for another generation is dropped. ``advance_generation`` marks a workspace
  boundary (a session load).
* ``ctx.restore_workspace`` - load a saved workspace slot as one document operation
  (installed by main; see ``session.py`` for the schema and the load policy).
* ``ctx.busy``    - raise/lower the footer progress bar without deferring anything, for
  work that already owns a tick of its own (``run_busy`` is this plus the deferral).
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
        self.submit_workflow = None  # (name, request, on_done) -> WorkflowHandle; by main
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
