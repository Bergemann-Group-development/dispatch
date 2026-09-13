"""Mutual-exclusion registry for the tools that own the right sidebar.

The range, normalise, filter, line-mask, atom and the three analysis tools (registered
by ``app.build``) all share the right sidebar and the canvas select gesture, so at most
one may be active at a time. ``ToolManager`` owns that invariant, so a tool never has to
know which others exist. A *tool* is any object exposing:

* ``active``      — bool, True while it owns the sidebar
* ``activate()``  — enter the tool (may no-op and leave ``active`` False if a
                    precondition fails, e.g. no spectrum selected)
* ``deactivate()``— leave the tool (idempotent)

The manager reads each tool's own ``active`` flag as the source of truth, so a
deferred activation (e.g. Filter's ``activate`` is scheduled on the next document
tick via ``run_busy`` so the progress bar paints first) stays correct: the flag is
only re-read on the next user action, by which time the deferred work has run.
"""
from __future__ import annotations


class ToolManager:
    def __init__(self):
        self._tools = {}     # name -> tool
        self._runners = {}   # name -> callable(activate_fn) wrapping the activation, if any

    def register(self, name, tool, runner=None):
        """Add a mutually-exclusive tool. ``runner`` (optional) wraps the activation,
        e.g. ``run_busy`` to paint the footer progress bar before a heavy enter."""
        self._tools[name] = tool
        if runner is not None:
            self._runners[name] = runner
        return tool

    def active_name(self):
        """The name of the tool that currently owns the sidebar, or ``None``."""
        return next((name for name, t in self._tools.items() if getattr(t, "active", False)),
                    None)

    def deactivate_active(self):
        """Stand down whichever tool is active (no-op if none is)."""
        name = self.active_name()
        if name is not None:
            self._tools[name].deactivate()

    def activate(self, name):
        """Activate ``name``, first standing down any other active tool."""
        if self.active_name() == name:
            return
        tool = self._tools[name]
        if name in ("filter", "normalise") and tool.canvas.selected() is None:
            tool.notify("Open and select a spectrum first.")
            return
        self.deactivate_active()
        runner = self._runners.get(name, lambda fn: fn())
        runner(self._tools[name].activate)

    def toggle(self, name):
        """Toggle ``name``: deactivate it if it is active, else activate it (standing
        down whatever else was active)."""
        if self.active_name() == name:
            self._tools[name].deactivate()
        else:
            self.activate(name)
