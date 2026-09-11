"""Registered atmosphere grids, edited in the Synthesis preferences page.

Each grid is ``{"id", "name", "folder", "pattern", "format"}``; a snapshot grid (a 3D
format) also carries ``snapshots`` and an optional ``resolution``, see
:mod:`tumult.synthesis.snapshot_grid`. The pattern is matched against each file's path
within the folder tree; ``tumult.synthesis.FORMAT_PATTERNS`` holds the stock pattern each
format starts from. The list
lives under the ``synthesis`` store section key ``atmos_grids`` (with ``atmos_grid_id``
naming the one selected in the Synth panel). This component renders an editable card per
grid — name / folder (+ an embedded :class:`FilePicker` to browse to it) / format / the
layout controls the format needs — plus an "add grid" button, and persists every edit.

Mirrors :class:`atompicker.AtomPicker`: a server-rendered list pushed through an
``html_bridge``, a hidden ``TextInput`` action channel for delegated clicks, a
``class_bridge`` to reveal the inline picker, and an embedded :class:`FilePicker`. Grid rows
are dynamic (ids vary), so the controls are raw ``<input>``/``<select>`` reporting through a
single ``grid-fields-input`` channel rather than the static :class:`controls.Form`.
"""
from __future__ import annotations

import html
import itertools
from pathlib import Path
import re

from bokeh.models import TextInput

import store
from tumult.synthesis import workflows as synth
from tumult.synthesis.snapshot_grid import DEFAULT_SNAPSHOTS, is_snapshot_grid
from filepicker import FilePicker
from htmlbridge import class_bridge, html_bridge, pulse

# atmosphere formats a grid may declare (the Fortran atmos_format values)
ATMOS_FORMATS = ["Marcs", "Stagger", "Stagger2", "Co5bold", "M3DIS", "Text", "save"]


def _new_id(existing):
    """A short stable id (``g1``, ``g2``, …) not already in ``existing``."""
    for n in itertools.count(1):
        cand = f"g{n}"
        if cand not in existing:
            return cand


class GridManager:
    """Owns the registered-grid list, its rendering, and the inline folder picker.

    ``on_change`` (optional) is called after any change so the panel can refresh its grid
    dropdown / file-existence warnings. ``start_dir`` is where the folder picker first opens.
    """

    def __init__(self, default_pattern, on_change=None, start_dir="~", builtin_dir=None):
        self.default_pattern = default_pattern
        self.on_change = on_change
        self.grids = [dict(g) for g in store.get_section("synthesis").get("atmos_grids", [])]
        self._register_builtin_marcs(builtin_dir)
        self._editing = None  # id of the grid whose folder the picker is choosing

        self.list_bridge = html_bridge("grid-list", name="grid-list-render")
        self.open_bridge = class_bridge("#grid-picker", "open")
        self.fields_input = TextInput(value="", visible=False, name="grid-fields-input")
        self.actions = TextInput(value="", visible=False, name="grid-actions")
        self.fields_input.on_change("value", self._on_field)
        self.actions.on_change("value", self._on_action)

        self.picker = FilePicker(self._on_folder_pick, container_id="grid-results",
                                 start_dir=start_dir)
        self.picker.layout.name = "grid_picker"

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [self.list_bridge, self.open_bridge, self.fields_input, self.actions]

    @property
    def results_name(self):
        return self.picker.results_name

    # --- built-in grids ---------------------------------------------------
    def _register_builtin_marcs(self, builtin_dir):
        """Auto-register the bundled MARCS grid the first time it appears on disk.

        Once the input data is fetched, ``input/atmos/MARCS`` holds the standard MARCS
        model-atmosphere grid. Rather than make every user hand-configure it in the Synthesis
        preferences, seed a ready-to-use grid entry pointing at that folder (using the MARCS
        filename ``default_pattern`` and ``Marcs`` format). Since it becomes the first
        registered grid, :meth:`synthsetup.SynthSetup.selected_grid_id` self-heals to select
        it, so synthesis works out of the box.

        Runs during construction, so it must NOT call ``on_change`` (the owning panel isn't
        built yet) — it writes the grid list to the store directly. A one-shot
        ``marcs_autoregistered`` flag records that we've seeded it, so a user who later deletes
        the grid does not get it silently re-created. The flag is only set once the folder
        actually exists, so fetching the data after first launch still triggers registration.
        """
        if not builtin_dir:
            return
        if store.get_section("synthesis").get("marcs_autoregistered"):
            return
        folder = (Path(builtin_dir).expanduser() / "MARCS")
        if not folder.is_dir():
            return  # input data not fetched yet; retry on a later launch
        store.update_pref("synthesis", "marcs_autoregistered", True)
        resolved = folder.resolve()
        already = any(g.get("id") == "marcs" for g in self.grids) or any(
            (g.get("folder") or "").strip()
            and Path(g["folder"]).expanduser().resolve() == resolved
            for g in self.grids)
        if already:
            return
        self.grids.append({"id": "marcs", "name": "MARCS", "folder": str(folder),
                           "pattern": self.default_pattern, "format": ATMOS_FORMATS[0]})
        store.update_pref("synthesis", "atmos_grids", [dict(g) for g in self.grids])

    # --- lookups ----------------------------------------------------------
    def _find(self, gid):
        return next((g for g in self.grids if g.get("id") == gid), None)

    @property
    def names(self):
        """``[(id, name)]`` for the panel's grid dropdown."""
        return [(g["id"], g.get("name") or g["id"]) for g in self.grids]

    @property
    def ref_grid_id(self):
        """Id of the grid flagged as the 1D LTE reference (for abundance corrections).

        Validated against the live registry; falls back to the first registered grid (and
        persists that fallback) so a reference is always defined while any grid exists."""
        ids = [g["id"] for g in self.grids]
        cur = store.get_section("synthesis").get("ref_grid_id")
        if cur not in ids:
            cur = ids[0] if ids else None
            if cur is not None:
                store.update_pref("synthesis", "ref_grid_id", cur)
        return cur

    def reference_grid(self):
        """The flagged 1D LTE reference grid dict, or ``None`` when no grids are registered."""
        return self._find(self.ref_grid_id)

    def missing_folders(self):
        """Names of grids whose configured folder does not exist (empty folders ignored)."""
        out = []
        for g in self.grids:
            folder = (g.get("folder") or "").strip()
            if folder and not Path(folder).expanduser().is_dir():
                out.append(g.get("name") or g["id"])
        return out

    # --- rendering --------------------------------------------------------
    def _card_html(self, g):
        gid = html.escape(g["id"])
        folder = g.get("folder") or ""
        warn = bool(folder.strip()) and not Path(folder).expanduser().is_dir()
        fmt = g.get("format") or ATMOS_FORMATS[0]
        opts = "".join(
            f'<option value="{html.escape(f)}"'
            f'{" selected" if f == fmt else ""}>{html.escape(f)}</option>'
            for f in ATMOS_FORMATS)
        folder_cls = "field-input grid-folder" + (" field-warn-input" if warn else "")
        return (
            f'<div class="grid-card" data-grid-id="{gid}">'
            f'<div class="grid-card-head">'
            f'<input type="text" class="field-input grid-name" data-grid-id="{gid}" '
            f'data-grid-key="name" value="{html.escape(g.get("name") or "")}" '
            f'placeholder="Grid name">'
            f'<button class="icon-btn grid-del" data-grid-id="{gid}" '
            f'title="Remove grid"><i class="ti ti-trash"></i></button></div>'
            f'<div class="path-row">'
            f'<input type="text" class="{folder_cls}" data-grid-id="{gid}" '
            f'data-grid-key="folder" value="{html.escape(folder)}" '
            f'placeholder="Atmosphere folder">'
            f'<button class="icon-btn grid-browse" data-grid-id="{gid}" '
            f'title="Browse for the grid folder"><i class="ti ti-folder-open"></i></button>'
            f'</div>'
            f'<select class="field-input grid-format" data-grid-id="{gid}" '
            f'data-grid-key="format">{opts}</select>'
            f'{self._layout_html(g, gid)}'
            f'</div>')

    @staticmethod
    def _layout_html(g, gid):
        """The controls below the format, which decide how the folder is read: the pattern
        binding stellar parameters (and, for snapshot grids, the snapshot number) from each
        file's path within the folder tree; for snapshot grids also which snapshots
        represent a node and the mesh their cubes are sampled onto."""
        out = [f'<input type="text" class="field-input grid-pattern" data-grid-id="{gid}" '
               f'data-grid-key="pattern" value="{html.escape(g.get("pattern") or "")}" '
               f'placeholder="Path pattern, e.g. {html.escape(synth.DEFAULT_PATTERN)}" '
               f'title="Matched against each file path relative to the folder; see '
               f'tumult.synthesis.compile_pattern for the token grammar">']
        if is_snapshot_grid(g):
            resolution = ",".join(str(v) for v in (g.get("resolution") or ()))
            out.append(
                f'<input type="text" class="field-input grid-snapshots" data-grid-id="{gid}" '
                f'data-grid-key="snapshots" '
                f'value="{html.escape(str(g.get("snapshots") or DEFAULT_SNAPSHOTS))}" '
                f'placeholder="Snapshots per node: a count, all, or numbers 8,31" '
                f'title="How many snapshots (in snapshot order), \'all\', or explicit '
                f'snapshot numbers represent each node; node spectra are their mean">'
                f'<input type="text" class="field-input grid-resolution" data-grid-id="{gid}" '
                f'data-grid-key="resolution" value="{html.escape(resolution)}" '
                f'placeholder="Mesh nx,ny,nz (blank: Engine preferences)">')
        return "".join(out)

    def _ref_html(self):
        """The '1D LTE reference grid' selector, shown above the cards when grids exist.

        Names the registered grid used as the LTE baseline by the NLTE-correction modal.
        Hidden entirely when no grids are registered (there is nothing to choose)."""
        if not self.grids:
            return ""
        ref = self.ref_grid_id
        opts = "".join(
            f'<option value="{html.escape(gid)}"{" selected" if gid == ref else ""}>'
            f'{html.escape(name)}</option>' for gid, name in self.names)
        return ('<div class="field-group"><label>1D LTE reference grid</label>'
                f'<select class="field-input grid-ref-select">{opts}</select></div>')

    def list_html(self):
        cards = "".join(self._card_html(g) for g in self.grids)
        if not cards:
            cards = '<div class="field-static">No grids registered.</div>'
        return (self._ref_html() + cards + '<button class="text-btn grid-add">'
                '<i class="ti ti-plus"></i> Add grid</button>')

    def render(self):
        self.list_bridge.value = self.list_html()

    # --- state ------------------------------------------------------------
    def _persist(self):
        store.update_pref("synthesis", "atmos_grids", [dict(g) for g in self.grids])
        if self.on_change:
            self.on_change()

    def _on_field(self, attr, old, new):
        if not new:
            return
        msg = new.rsplit(" ", 1)[0]              # strip the trailing nonce
        ident, _, value = msg.partition("=")     # value may contain spaces / "{...}"
        gid, _, key = ident.partition("|")
        g = self._find(gid)
        if g is None or key not in ("name", "folder", "pattern", "format", "snapshots",
                                    "resolution"):
            return
        if key == "resolution":
            parts = [p for p in re.split(r"[,\sx]+", value.strip()) if p]
            if parts and not (len(parts) == 3 and all(p.isdigit() for p in parts)):
                self.render()       # reject silently by restoring the stored value
                return
            g[key] = [int(p) for p in parts] if parts else None
        else:
            g[key] = value
            if key == "format" and (g.get("pattern") or "") in ("", *synth.FORMAT_PATTERNS.values()):
                # a format switch reseeds a pattern the user has not customised with the
                # format's stock pattern (tumult.synthesis.FORMAT_PATTERNS)
                g["pattern"] = synth.default_pattern(value)
        self._persist()
        # name/folder edits change the dropdown label / warning state, and the format
        # decides which layout controls a card shows -> re-render the cards
        if key in ("name", "folder", "format"):
            self.render()

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]             # strip the trailing nonce
        if kind == "add":
            gid = _new_id({g["id"] for g in self.grids})
            self.grids.append({"id": gid, "name": "New grid", "folder": "",
                               "pattern": self.default_pattern, "format": ATMOS_FORMATS[0]})
            self._persist()
            self.render()
        elif kind == "del":
            self.grids = [g for g in self.grids if g.get("id") != arg]
            self._persist()
            self.render()  # ref_grid_id self-heals to a live grid via _ref_html on render
        elif kind == "setref":
            if self._find(arg) is not None:
                store.update_pref("synthesis", "ref_grid_id", arg)
                if self.on_change:
                    self.on_change()
        elif kind == "browse":
            self._editing = arg
            g = self._find(arg)
            start = (g.get("folder") or "").strip() if g else ""
            self.picker.query.value = ""
            self.picker._scan(f"{Path(start).expanduser()}/" if start else "")
            pulse(self.open_bridge, "on")
        elif kind == "folder":
            # "Select this folder": adopt the picker's current directory for the editing grid
            self._adopt_folder(self.picker.current_dir)

    def _on_folder_pick(self, filepath):
        # FilePicker only fires on_select for files; folders navigate. Allow picking a file
        # too — adopt its parent directory as the grid folder.
        self._adopt_folder(Path(filepath).parent)

    def _adopt_folder(self, folder):
        g = self._find(self._editing)
        if g is None or folder is None:
            return
        g["folder"] = str(folder)
        self._persist()
        self.render()
        pulse(self.open_bridge, "off")
