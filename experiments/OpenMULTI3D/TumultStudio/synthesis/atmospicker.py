"""The atmosphere picker of a synthesis setup: a :class:`filepicker.FilePicker` that lists
the parsed models of a registered grid instead of the files of its folder.

A grid's models are what its pattern matched (:func:`tumult.synthesis.discover_atmos_grid`,
the same nodes interpolation uses), one row per model or, for a snapshot grid, per
snapshot. A row reads ``5750 K · log g 4.50 · [Fe/H] +0.00 · ξ 1`` and the query box is a
filter over that text: every space-separated fragment has to occur, so ``5750 4.5`` or
``-2.00 giant`` narrows a fifteen-thousand-model MARCS grid to a handful. The chosen row's
label then shows in the box, as a picked path does in file mode.

``set_root`` puts the picker back in file mode (the Individual file), listing a
directory as :class:`FilePicker` does. A grid is scanned away from the document, and the
picker is in node mode from the moment the scan starts, so the box says what it is for
while a large grid is read.
"""
from __future__ import annotations

import html
from pathlib import Path

from tumult.synthesis import workflows as synth
from tumult.synthesis.snapshot_grid import is_snapshot_grid
from ..ui.filepicker import FilePicker
from ..ui.latestwork import LatestWork

FILE_PLACEHOLDER = "Type a path…"
NODE_PLACEHOLDER = "Filter models, e.g. 5750 4.5 -1"
SCANNING = "Listing the grid's models…"


def node_label(teff, logg, feh, vmic=None, snap=None):
    """The row text of one model."""
    out = f"{teff:.0f} K · log g {logg:.2f} · [Fe/H] {feh:+.2f}"
    if vmic is not None:
        out += f" · ξ {vmic:g}"
    if snap is not None:
        out += f" · snapshot {int(snap)}"
    return out


def grid_nodes(grid, root=None):
    """``[(label, path, sort key)]`` for every model of a registered grid record, every
    snapshot of a snapshot grid listed on its own. Raises ``FileNotFoundError`` for a
    missing folder and ``ValueError`` when nothing matched the pattern."""
    folder = Path(grid["folder"]).expanduser()
    if root is not None and not folder.is_absolute():
        folder = Path(root) / folder
    pattern = grid.get("pattern") or synth.DEFAULT_PATTERN
    compiled = synth.compile_pattern(pattern)
    snapshots = "all" if is_snapshot_grid(grid) else None
    nodes = synth.discover_atmos_grid(folder, pattern, snapshots=snapshots)
    rows = []
    for key, value in nodes.items():
        teff, logg, feh = key[:3]
        vmic = key[3] if len(key) > 3 else None
        for path in (value if isinstance(value, (list, tuple)) else (value,)):
            snap = None
            if snapshots is not None:
                found = compiled.search(Path(path).relative_to(folder).as_posix())
                snap = found.get("snap") if found else None
            rows.append((node_label(teff, logg, feh, vmic, snap), path,
                         (teff, logg, feh, vmic or 0.0, snap or 0, Path(path).name)))
    rows.sort(key=lambda r: r[2])
    return rows


class AtmosPicker(FilePicker):
    def __init__(self, on_select, container_id, root=None, max_items=300):
        self._nodes = None        # [(label, path, sort key)] in node mode, else None
        self._problem = None      # why the node list is empty, when it is
        self._by_path = {}
        self._note = ""
        super().__init__(on_select, container_id, root=root, max_items=max_items)
        self._grid_work = LatestWork(self._submit_scan, "discover atmosphere grid")

    # --- modes ----------------------------------------------------------------
    def set_root(self, root):
        """File mode: list ``root`` (the Individual file folder)."""
        if hasattr(self, "_grid_work"):
            self._grid_work.invalidate()
        self._nodes, self._problem, self._by_path = None, None, {}
        self.query.placeholder = FILE_PLACEHOLDER
        super().set_root(root)

    def set_nodes(self, rows, problem=None):
        """Node mode: filter ``rows`` (see :func:`grid_nodes`); ``problem`` is shown in place
        of an empty list."""
        self._nodes = list(rows)
        self._problem = problem
        self._by_path = {self._key(path): label for label, path, _ in self._nodes}
        self.query.placeholder = NODE_PLACEHOLDER
        self.root = None
        self.query.value = ""
        self._scan("")

    def set_grid(self, grid, root=None, on_ready=None):
        """Discover a grid away from the document and install only its latest result.

        The picker is put in node mode at once, listing nothing but saying the scan is
        under way, so a visible picker never shows the file browser of another mode while
        a grid is being read. ``on_ready`` is called on the document thread once the models
        are installed (the caller's chance to show the picked one by its label, which only
        exists once they are listed).
        """
        marker = (dict(grid), root)
        if not self.node_mode:
            self.set_nodes([], SCANNING)

        def compute(token):
            try:
                return grid_nodes(marker[0], marker[1]), None
            except (OSError, ValueError) as exc:
                return [], f"Cannot list the grid: {exc}"

        def apply(result):
            rows, problem = result
            self.set_nodes(rows, problem)
            if on_ready is not None:
                on_ready()

        self._grid_work.request(compute, apply)

    @property
    def node_mode(self):
        return self._nodes is not None

    def label_for(self, path):
        """The row label of a model path, or ``None`` outside node mode / off the list."""
        return self._by_path.get(self._key(path)) if path else None

    # --- listing ----------------------------------------------------------------
    def _scan(self, text):
        if self._nodes is None:
            self._note = ""
            super()._scan(text)
            return
        self.current_dir = None
        if self._problem:
            self._note = ""
            self._list([], self._problem)
            return
        wanted = [w.lower() for w in (text or "").split()]
        hits = [(label, path) for label, path, _ in self._nodes
                if all(w in label.lower() for w in wanted)]
        shown = hits[:self.max_items]
        self._note = (f"Showing {len(shown)} of {len(hits)} models; refine the filter"
                      if len(hits) > len(shown) else "")
        self._list([(label, path, False, "file-text") for label, path in shown],
                   "No model matches the filter")

    def _render(self):
        super()._render()
        if self._note:
            self.results.value += f'<div class="file-empty">{html.escape(self._note)}</div>'

    def _on_pick(self, attr, old, new):
        if self._nodes is None or not new:
            super()._on_pick(attr, old, new)
            return
        self._load(Path(new.split("\x00", 1)[0]))
