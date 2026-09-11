"""Snapshot grids: 3D model grids whose nodes are time series of snapshots.

Discovery is the same recursive pattern scan as for a 1D grid
(:func:`~tumult.synthesis.workflows.discover_atmos_grid`); what differs is that a node
maps to several files, and which of them represent the node is a user setting stored on
the registered grid record as ``snapshots``. The Fortran reads one snapshot per run, so a
node's spectrum is the mean over its selected snapshots (see ``api._collapse_snapshots``),
and the node's velocity field replaces microturbulence, so snapshot nodes carry no vmic.
"""
from __future__ import annotations

from pathlib import Path
import re

# atmos_format values (lower case, as the Fortran compares) whose grids hold snapshots
SNAPSHOT_FORMATS = ("stagger", "stagger2", "co5bold", "m3dis", "must", "multi")
DEFAULT_SNAPSHOTS = "1"


def parse_snapshot_selection(value):
    """A snapshot policy: ``'all'``, a count ``'3'`` (the first n in snapshot order) or
    explicit snapshot numbers ``'8,31,204'``. Returns ``('all', None)``, ``('count', n)``
    or ``('numbers', {int, ...})``."""
    text = str(value if value is not None else DEFAULT_SNAPSHOTS).strip().lower()
    text = text or DEFAULT_SNAPSHOTS
    if text == "all":
        return "all", None
    parts = [p.strip() for p in re.split(r"[,\s]+", text) if p.strip()]
    if len(parts) == 1 and parts[0].isdigit():
        count = int(parts[0])
        if count < 1:
            raise ValueError("Snapshot count must be at least 1.")
        return "count", count
    numbers = set()
    for part in parts:
        if not part.lstrip("#").isdigit():
            raise ValueError(f"Snapshot selection {value!r} is not 'all', a count or "
                             "comma-separated snapshot numbers.")
        numbers.add(int(part.lstrip("#")))
    return "numbers", numbers


def is_snapshot_grid(grid):
    """Whether a registered grid record describes a grid of 3D snapshots."""
    return str((grid or {}).get("format") or "").lower() in SNAPSHOT_FORMATS


def is_snapshot_node(atmos_files):
    """Whether a discovered grid maps nodes to snapshot sequences rather than files."""
    return any(isinstance(value, (list, tuple)) for value in atmos_files.values())


def discover_grid(grid, root=None):
    """Discover a registered grid record: ``folder`` and ``pattern``, plus the
    ``snapshots`` policy when its ``format`` is a snapshot format. ``root`` resolves a
    relative folder."""
    from .workflows import DEFAULT_PATTERN, discover_atmos_grid
    folder = Path(grid["folder"]).expanduser()
    if root is not None and not folder.is_absolute():
        folder = Path(root) / folder
    pattern = grid.get("pattern") or DEFAULT_PATTERN
    if is_snapshot_grid(grid):
        return discover_atmos_grid(folder, pattern,
                                   snapshots=grid.get("snapshots") or DEFAULT_SNAPSHOTS)
    return discover_atmos_grid(folder, pattern)
