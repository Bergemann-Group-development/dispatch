"""Snapshot grids: 3D model grids whose nodes are time series of snapshots.

Discovery is the same recursive pattern scan as for a 1D grid
(:func:`~tumult.synthesis.workflows.discover_atmos_grid`); what differs is that a node
maps to several files, and which of them represent the node is a user setting stored on
the registered grid record as ``snapshots``. The Fortran reads one snapshot per run, so a
node's spectrum is the mean over its selected snapshots (see ``api._collapse_snapshots``),
and the node's velocity field replaces microturbulence, so snapshot nodes carry no vmic.

A count policy picks snapshots in :func:`halton_order`, a fixed Owen-scrambled base-2
Halton ordering of the node's available snapshots. Its prefixes nest (the ``k`` picks are
among the ``k + 1`` picks) and spread over the time series without the equal spacing that
would alias against oscillations, and a grid download with ``k`` snapshots per model
(:mod:`.grid_catalog`) fetches exactly the snapshots ``snapshots="k"`` then selects.
"""
from __future__ import annotations

from pathlib import Path
import re

# atmos_format values (lower case, as the Fortran compares) whose grids hold snapshots
SNAPSHOT_FORMATS = ("stagger", "stagger2", "co5bold", "m3dis", "must", "multi")
DEFAULT_SNAPSHOTS = "1"


# the seed of the Owen scramble; fixed so that a download and a later selection agree
HALTON_SEED = 0x9E3779B9


def _reverse_bits(x):
    return int(f"{x & 0xFFFFFFFF:032b}"[::-1], 2)


def _laine_karras(x, seed):
    """Laine-Karras hash permutation of a 32-bit word (Burley 2020, "Practical
    hash-based Owen scrambling"): every output bit depends only on its own and lower
    input bits, which on a bit-reversed radical inverse is exactly the nesting an Owen
    scramble needs."""
    x = (x + seed) & 0xFFFFFFFF
    for c in (0x6C50B47C, 0xB82F1E52, 0xC7AFE638, 0x8D22F6E6):
        x ^= (x * c) & 0xFFFFFFFF
    return x & 0xFFFFFFFF


def scrambled_halton(index, seed=HALTON_SEED):
    """Term ``index`` (0-based) of the Owen-scrambled base-2 Halton sequence, in [0, 1).

    The unscrambled sequence is the radical inverse of ``index``; the scramble permutes
    every dyadic interval's halves by a hash of the digits before it. Any ``2**m`` leading
    terms still hold one point per interval of width ``2**-m``, but the points sit at
    hashed positions within their intervals rather than at equal spacing."""
    return _reverse_bits(_laine_karras(index & 0xFFFFFFFF, seed)) / 2.0 ** 32


def halton_order(n, seed=HALTON_SEED):
    """The indices ``0 .. n-1`` in the order the scrambled Halton sequence visits them.

    Each term maps to ``floor(u * n)``; an index already visited is skipped, so the
    result is a permutation of ``range(n)`` whose every prefix is a well-spread subset.
    Selecting ``k`` snapshots means taking the first ``k`` entries."""
    order = []
    seen = set()
    index = 0
    while len(order) < n:
        i = int(scrambled_halton(index, seed) * n)
        if i not in seen:
            seen.add(i)
            order.append(i)
        index += 1
    return order


def select_count(available, count, seed=HALTON_SEED):
    """The first ``count`` of ``available`` (a sequence in snapshot order) in
    :func:`halton_order`, returned in snapshot order; all of them when ``count``
    exceeds the sequence."""
    items = list(available)
    picks = sorted(halton_order(len(items), seed)[:count])
    return [items[i] for i in picks]


def parse_snapshot_selection(value):
    """A snapshot policy: ``'all'``, a count ``'3'`` (the first n in :func:`halton_order`) or
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
