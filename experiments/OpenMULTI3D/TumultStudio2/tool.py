"""The base class of the tools that take turns owning the right sidebar.

The canvas gestures (``RangeTool`` for the Select range tool's range, ``FilterTool``,
``NormaliseTool``), the object views (``MaskTool``, ``AtomTool``) and the three
analysis wizards (``AbundanceTool``, ``AbundCorrTool``, ``IonBalanceTool``) share one
protocol: ``active``, ``activate()`` and ``deactivate()``, which ``ToolManager`` drives
for mutual exclusion; a panel bridge for the sidebar body; on/off class bridges for the
sidebar mode and, where the tool has one, its toolbar button; and ``update()``, called
on the context's ``"changed"`` event when the workspace changes. ``Tool`` owns that
plumbing. Each tool keeps the activation logic and is the source of truth for its own
``active`` flag, including a deferred entry (Filter activates on the next tick) and a
refused one (a missing precondition leaves the flag down).

A tool that edits workspace state in place (``NormaliseTool``) runs an *edit session*:
it calls :meth:`begin_edit` on entry, which keeps what :meth:`capture_edit` returns.
Leaving the tool by any route keeps the edits (Done, its toolbar button, switching tools,
selecting another object); :meth:`cancel` (the panel's Cancel and Esc, which the manager
routes to tools with ``cancellable`` set) hands that snapshot to :meth:`restore_edit`
before leaving.

Subclasses set ``name`` and call :meth:`own_mode_bridges` / :meth:`own_panel` in
``__init__``; ``push_mode`` / ``push_panel`` then drive them. :func:`suppress_callbacks`
is the shared guard for setters whose Bokeh writes would otherwise re-enter a drag handler.
"""
from __future__ import annotations

from functools import wraps

from htmlbridge import pulse


def suppress_callbacks(method):
    """Run a synchronous setter with ``self._suppress`` raised, so the Bokeh change
    callbacks its writes trigger (a knot or crosshair drag handler) do not refit.

    The previous flag is restored in ``finally``, so a failing write cannot leave later
    user drags ignored and a nested suppressed call leaves the outer one suppressed."""
    @wraps(method)
    def wrapped(self, *args, **kwargs):
        previous = self._suppress
        self._suppress = True
        try:
            return method(self, *args, **kwargs)
        finally:
            self._suppress = previous
    return wrapped


class Tool:
    name = ""
    cancellable = False     # Esc reaches cancel() (see ToolManager.cancel_active)

    def __init__(self, ctx):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.active = False
        self._mode_bridges = ()
        self._panel_bridge = None
        self._edit_baseline = None
        ctx.subscribe("changed", lambda **_: self.update())

    # notify always resolves through the context, so tools built before the footer
    # exists still pick up the real notifier once main installs it.
    @property
    def notify(self):
        return self.ctx.notify

    # --- bridge ownership ----------------------------------------------------
    def own_mode_bridges(self, *pairs):
        """Create the on/off class bridges this tool pulses (e.g. its right-sidebar
        panel visibility and its toolbar button state). ``pairs`` are ``(selector,
        klass)`` tuples."""
        self._mode_bridges = tuple(self.ctx.bridges.cls(sel, klass) for sel, klass in pairs)

    def own_panel(self, container_id):
        """Create the html bridge for this tool's inspector-panel body."""
        self._panel_bridge = self.ctx.bridges.html(container_id)

    def push_mode(self, state):
        """Pulse "on"/"off" to all of the tool's mode bridges."""
        for bridge in self._mode_bridges:
            pulse(bridge, state)

    def push_panel(self):
        """Re-render the tool's panel HTML into its panel bridge."""
        if self._panel_bridge is not None:
            self._panel_bridge.value = self.panel_html()

    # --- protocol (subclass responsibilities) ---------------------------------
    def panel_html(self):
        """The inspector-panel body ("" when inactive)."""
        return ""

    def activate(self):
        """Enter the tool; may refuse (leave ``active`` False) on a failed
        precondition. Subclasses implement."""
        raise NotImplementedError

    def deactivate(self):
        """Leave the tool (idempotent). Subclasses implement."""
        raise NotImplementedError

    # --- edit sessions ----------------------------------------------------------
    def begin_edit(self):
        """Snapshot the state the tool is about to edit."""
        self._edit_baseline = self.capture_edit()

    def end_edit(self):
        """Forget the snapshot: the edits are kept."""
        self._edit_baseline = None

    def cancel(self):
        """Leave the tool with the edited state as it was on entry."""
        baseline, self._edit_baseline = self._edit_baseline, None
        if self.active and baseline is not None:
            self.restore_edit(baseline)
        self.deactivate()

    def capture_edit(self):
        """The state an edit session may change, as :meth:`restore_edit` takes it."""
        return None

    def restore_edit(self, baseline):
        """Put back the state :meth:`capture_edit` returned."""

    def update(self):
        """The workspace changed (``"changed"`` event); resync or stand down."""
