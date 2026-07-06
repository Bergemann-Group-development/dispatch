"""The formal base class for the mutually-exclusive sidebar tools.

The EW / Filter / Normalise / Synth / Abundance tools always shared an implicit
protocol - ``active``, ``activate()``, ``deactivate()`` (what ``ToolManager``
duck-types), a panel bridge, a pair of on/off mode bridges, and an ``update()``
called when the workspace changes. ``Tool`` makes that protocol explicit and owns the
common plumbing; the activation/deactivation *logic* stays with each tool (it is the
source of truth for its own ``active`` flag, including deferred or refused entries).

Subclasses set ``name`` and typically call :meth:`own_mode_bridges` /
:meth:`own_panel` in ``__init__``; ``push_mode``/``push_panel`` then drive them.
``update()`` is subscribed to the context's ``"changed"`` event automatically.
"""
from __future__ import annotations

from htmlbridge import pulse


class Tool:
    name = ""

    def __init__(self, ctx):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.active = False
        self._mode_bridges = ()
        self._panel_bridge = None
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

    def update(self):
        """The workspace changed (``"changed"`` event); resync or stand down."""
