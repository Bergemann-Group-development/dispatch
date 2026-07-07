"""Shared building blocks for the sortable HTML tables.

The Line Navigator (:mod:`linenav`) and the two abundance wizards (:mod:`abundancetool`
/ :mod:`abundcorrtool`) all render a ``<table class="linelist-table">`` whose column
headers are clickable to sort. This module owns the sort *state* and the header markup so
those tools stay in sync; each tool still supplies its own row rendering and its own data
ordering (a DataFrame for the line list, a list of line objects for the abundance tables).

The client echoes a clicked header back through each tool's own action channel as the
column key in the ``data-sortcol`` attribute; the tool feeds it to :meth:`SortState.toggle`
and re-renders.
"""
from __future__ import annotations

import html
import math


class SortState:
    """Which column a table is sorted by, and in which direction.

    ``columns`` maps each sortable column key to its default *ascending?* flag: clicking a
    fresh column adopts that default, clicking the already-active column flips direction
    (this matches how spreadsheet-style tables behave).
    """

    def __init__(self, columns, col, asc=None):
        self.columns = dict(columns)
        self.col = col
        self.asc = self.columns[col] if asc is None else asc

    def toggle(self, col):
        """Apply a header click. Returns True when the sort changed (caller re-renders)."""
        if col not in self.columns:
            return False
        if col == self.col:
            self.asc = not self.asc
        else:
            self.col, self.asc = col, self.columns[col]
        return True

    def order(self, items, key):
        """Indices of ``items`` in the active order; NaN/None keys always land last.

        ``key(item)`` yields the sort value (a number or a string). Numeric columns keep
        missing values (None / NaN) at the bottom regardless of direction — mirroring
        pandas' ``na_position="last"`` used by the DataFrame-backed line list.
        """
        idx = list(range(len(items)))
        vals = [key(items[i]) for i in idx]
        if vals and all(isinstance(v, str) for v in vals):
            idx.sort(key=lambda i: vals[i].lower(), reverse=not self.asc)
            return idx
        live = [i for i in idx if _present(vals[i])]
        dead = [i for i in idx if not _present(vals[i])]
        live.sort(key=lambda i: vals[i], reverse=not self.asc)
        return live + dead


def _present(v):
    """True unless ``v`` is a missing numeric value (None or NaN)."""
    if v is None:
        return False
    return not (isinstance(v, float) and math.isnan(v))


def sort_header(label, col, state):
    """A ``<th>`` for one column.

    ``label`` is trusted HTML (callers pass ``<sub>``/``<sup>`` markup). When ``col`` is a
    sortable key it becomes a clickable ``.line-sort`` header carrying ``data-sortcol`` and
    a direction arrow while active; a ``col`` of ``None`` (or one absent from ``state``)
    renders a plain, non-sortable header — used for the mask checkbox and flags columns.
    """
    if col is None or col not in state.columns:
        return f"<th>{label}</th>"
    col_attr = html.escape(str(col), quote=True)
    if col == state.col:
        arrow = "▲" if state.asc else "▼"
        return (f'<th class="line-sort active" data-sortcol="{col_attr}">'
                f"{label} {arrow}</th>")
    return f'<th class="line-sort" data-sortcol="{col_attr}">{label}</th>'


def sort_headers(state, columns):
    """Join :func:`sort_header` over a sequence of ``(label, col)`` pairs."""
    return "".join(sort_header(label, col, state) for label, col in columns)
