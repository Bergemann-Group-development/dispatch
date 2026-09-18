"""The local continuum of a single measured line, edited on the spectrum.

A measured equivalent width stands on a continuum. The global one the spectrum carries is
usually it, but a line sitting in a blend or on a badly placed knot needs its own: two nodes
at the edges of its synthesis window, which :func:`tumult.abundance.set_local_continuum`
re-measures the line through. That per-line choice is recorded on the line (``renorm`` and
``continuum``), so it travels with the result, and it is what the ``Local continuum``
checkbox column switches.

:class:`LocalContinuum` owns the interaction: the two draggable markers and the dashed line
between them on the main figure, the debounce that measures only once a drag has settled,
and the delivery of the re-measured line back into the result. Every measurement runs off the
document thread through :class:`~..ui.latestwork.LatestWork` on a deep copy, so a
failure leaves the displayed result exactly as it was and only the latest request is applied.

The view that owns the result supplies what this cannot know: which spectrum the
result was measured on (``spec_of``), how its equivalent-width hatch is drawn through the
edited continuum (``draw_hatch``), whether a result is still the listed object of its uid
once a measurement settles (``listed``), which line is in view (``current``) and what to
repaint when an edit lands (``on_edited``).
"""
from __future__ import annotations

from copy import deepcopy

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, PointDrawTool

from tumult import abundance as ab

from ..ui.latestwork import LatestWork
from ..ui.style import ACCENT
from ..ui.debounce import Trailing
from ..ui.units import to_value

# The "Local continuum" header of the per-line table: a checkbox column, see check_cols.
LOCAL_CONTINUUM_HEAD = "Local continuum"

DRAG_SETTLE_MS = 200         # quiet interval after a marker drag before it is measured


def continuum_cell(ln, index):
    """The per-line checkbox cell of the ``Local continuum`` column."""
    checked = " checked" if ln.renorm else ""
    return (f'<td><input type="checkbox" class="line-continuum-chk" '
            f'data-line="{index}"{checked} title="Adjust this measurement’s local continuum">'
            '</td>')


class LocalContinuum:
    """The local-continuum markers of the line currently in view, and their measurement."""

    def __init__(self, ctx, prefix, *, spec_of, draw_hatch, on_edited, current, listed):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self.spec_of = spec_of            # (result) -> Spectrum | None
        self.draw_hatch = draw_hatch      # (result, line) -> None: the EW hatch through it
        self.on_edited = on_edited        # (result, [line index]) -> None: repaint
        self.current = current            # () -> (result | None, jumped index | None)
        self.listed = listed              # (result) -> bool: still the object of its uid
        self._sync = False
        self.source = ColumnDataSource(dict(x=[], y=[]))
        # These are editable markers, not a spectrum: canvas.py keeps them off the main
        # figure's y auto-fit (see its note on ``main.y_range.renderers``), so dragging one
        # cannot widen the range to chase the very point being moved.
        self.markers = self.canvas.main.scatter(
            "x", "y", source=self.source, size=10, color=ACCENT, visible=False)
        self.curve = self.canvas.main.line(
            "x", "y", source=self.source, color=ACCENT, line_dash="dashed", visible=False)
        self.drag = PointDrawTool(renderers=[self.markers], add=False)
        self.canvas.main.add_tools(self.drag)
        self.source.on_change("data", self._on_drag)
        self._delay = Trailing(DRAG_SETTLE_MS)
        self._immediate = getattr(ctx.doc, "session_context", None) is None
        self._work = LatestWork(
            lambda *args: ctx.submit(*args), f"{prefix} continuum measurement",
            lambda exc: ctx.notify(f"Local continuum could not be measured: {exc}"),
            immediate=self._immediate)

    # --- the markers ---------------------------------------------------------------
    def hide(self):
        self.markers.visible = self.curve.visible = False
        self.canvas.release_edit_tool(self.drag)

    def render(self, res, index):
        """Show the nodes of line ``index`` of ``res``, or hide them when that line has no
        local continuum of its own."""
        ln = res.lines[index]
        if not ln.renorm or not ln.continuum:
            return self.hide()
        xx, yy = np.asarray(ln.continuum).T
        self._sync = True
        try:
            self.source.data = dict(x=to_value(xx, u.AA, self.canvas.xunit), y=yy)
        finally:
            self._sync = False
        self.markers.visible = self.curve.visible = True
        self.canvas.claim_edit_tool(self.drag)
        self.draw_hatch(res, ln)

    # --- editing -------------------------------------------------------------------
    def edit(self, res, index, **kwargs):
        """Re-measure line ``index`` of ``res`` under ``kwargs`` (``enabled=`` or
        ``nodes=``, as :func:`tumult.abundance.set_local_continuum` takes them)."""
        spec = self.spec_of(res)
        if spec is None:
            self.ctx.notify("The target spectrum is no longer loaded.")
            return

        def compute(token):
            token.check()
            working = deepcopy(res)
            ab.set_local_continuum(working, working.lines[index], spec, **kwargs)
            token.check()
            state = {key: value for key, value in working.__dict__.items() if key != "lines"}
            return working.lines[index], state

        def apply(result):
            if not self.listed(res) or index >= len(res.lines):
                return
            edited, state = result
            res.lines[index].__dict__.update(edited.__dict__)
            res.__dict__.update(state)
            self.on_edited(res, [index])

        self._work.request(compute, apply)

    def set_enabled(self, res, lines, enabled):
        """Switch the local continuum of several lines at once; each is remeasured."""
        spec = self.spec_of(res)
        if spec is None:
            self.ctx.notify("The target spectrum is no longer loaded.")
            return
        indices = [res.lines.index(ln) for ln in lines]

        def compute(token):
            token.check()
            working = deepcopy(res)
            edits = [(working.lines[index], {"enabled": enabled}) for index in indices]
            failed = ab.set_local_continua(working, edits, spec)
            failed_by_id = {id(line): exc for line, exc in failed}
            outcome = [(index, working.lines[index], failed_by_id.get(id(working.lines[index])))
                       for index in indices]
            token.check()
            state = {key: value for key, value in working.__dict__.items() if key != "lines"}
            return outcome, state

        def apply(result):
            if not self.listed(res):
                return
            outcome, state = result
            failures = []
            for index, edited, exc in outcome:
                if exc is not None:
                    failures.append(f"{res.lines[index].lam0:.3f} Å: {exc}")
                else:
                    res.lines[index].__dict__.update(edited.__dict__)
            res.__dict__.update(state)
            if failures:
                self.ctx.notify("Local continuum unchanged for " + "; ".join(failures))
            self.on_edited(res, indices)

        self._work.request(compute, apply)

    # --- dragging ------------------------------------------------------------------
    def _on_drag(self, attr, old, new):
        if self._sync or not self.markers.visible:
            return
        if self._immediate:
            self._commit_drag()
        else:
            self._delay("continuum-drag", self._commit_drag)

    def _commit_drag(self):
        """Measure only after the marker drag has been quiet for a short interval."""
        res, index = self.current()
        if res is None or index is None:
            return
        ln = res.lines[index]
        new = self.source.data
        if len(new.get("y", [])) != 2:
            self.render(res, index)
            return
        # The markers move vertically; their wavelengths remain at the synth edges.
        nodes = tuple(zip(ln.window, (float(y) for y in new["y"])))
        self.edit(res, index, nodes=nodes)
