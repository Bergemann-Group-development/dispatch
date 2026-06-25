"""A width-only box-select gesture over the canvas main figure, with a shaded band.

Shared by the EW and Synthesis tools. Each instance owns its own ``BoxSelectTool`` and
``BoxAnnotation``; ``_on_geometry`` only reacts when *its* tool is the active drag, so
several range-select tools can coexist on one figure without cross-talk. Subclasses
override ``on_select(x0, x1)`` (the band is shown for them by default).
"""
from __future__ import annotations

from bokeh.events import SelectionGeometry
from bokeh.models import BoxAnnotation, BoxSelectTool


class RangeSelect:
    def __init__(self, canvas, color, description, fill_alpha=0.07):
        self.canvas = canvas
        self.main = canvas.main
        self.selection = None  # (x0, x1) in display units, or None
        self.sel_box = BoxAnnotation(fill_color=color, fill_alpha=fill_alpha,
                                     line_color=color, line_alpha=0.4, visible=False)
        self.main.add_layout(self.sel_box)
        self.tool = BoxSelectTool(dimensions="width", description=description,
                                  continuous=False)
        self.main.add_tools(self.tool)
        self.main.on_event(SelectionGeometry, self._on_geometry)

    # --- tool arm / band --------------------------------------------------
    def arm(self):
        """Make this tool the figure's active drag gesture."""
        self.main.toolbar.active_drag = self.tool

    def disarm(self):
        """Release the drag gesture (back to panning) if we still hold it."""
        if self.main.toolbar.active_drag is self.tool:
            self.main.toolbar.active_drag = "auto"

    def show_band(self, x0, x1):
        """Record the drawn range and shade it; drop any stale point selection."""
        self.selection = (min(x0, x1), max(x0, x1))
        self.sel_box.left, self.sel_box.right = self.selection
        self.sel_box.visible = True
        for o in self.canvas.objects:
            o.main_source.selected.indices = []

    def clear_selection(self):
        self.selection = None
        self.sel_box.visible = False

    # --- gesture ----------------------------------------------------------
    def _on_geometry(self, event):
        if self.main.toolbar.active_drag is not self.tool:
            return  # another tool's gesture — ignore
        if not getattr(event, "final", True):
            return
        g = event.geometry
        if g.get("x0") is None or g.get("x1") is None:
            return
        self.on_select(float(g["x0"]), float(g["x1"]))

    def on_select(self, x0, x1):
        """Hook: a range was drawn (display units). Default just shades it."""
        self.show_band(x0, x1)
