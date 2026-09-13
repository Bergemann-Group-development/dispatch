"""Atmosphere grids, the Preferences page of that name.

The page holds one storage folder and one snapshot count for every grid, then the
published grids of :mod:`tumult.synthesis.grid_catalog`, one row each with what is on disk
and *Download all* / *Edit* / *Delete*; *Edit* lists the grid's models with a download or
delete button per model. Custom grids (a folder we do not host) keep an editable card:
name / folder (its field hands over to an embedded :class:`FilePicker`) / format / the
layout controls the format needs; ``tumult.synthesis.FORMAT_PATTERNS`` holds the stock
pattern each format starts from.

The registered grids other tools read (:attr:`grids`, :attr:`names`,
:meth:`reference_grid`) are the catalogue grids with models under the storage folder,
each as the record :func:`grid_catalog.grid_record` builds with the shared snapshot count,
plus the custom records. The store's ``synthesis`` section keeps ``grid_storage``,
``grid_snapshots``, the custom records under ``atmos_grids``, and ``atmos_grid_id`` /
``ref_grid_id`` naming the grid selected in the Synth panel and the 1D LTE reference.

Listing a grid and every download run as document jobs on the IO pool
(``ctx.submit_io``); progress lands in the page's status line. Mirrors
:class:`atompicker.AtomPicker`: a server-rendered page pushed through an ``html_bridge``,
a hidden ``TextInput`` action channel for delegated clicks, a ``class_bridge`` to reveal
the inline picker, and an embedded :class:`FilePicker`. Rows are dynamic, so the controls
are raw ``<input>``/``<select>`` reporting through the single ``grid-fields-input`` channel.
"""
from __future__ import annotations

import html
import itertools
from pathlib import Path
import re
import time

from bokeh.models import TextInput

import store
from tumult.synthesis import grid_catalog
from tumult.synthesis.cancel import CancelledError
from tumult.synthesis import workflows as synth
from tumult.synthesis.snapshot_grid import DEFAULT_SNAPSHOTS, is_snapshot_grid
from filepicker import FilePicker
from htmlbridge import class_bridge, html_bridge, pulse
from tableui import SortState, sort_headers

# atmosphere formats a custom grid may declare (the Fortran atmos_format values)
ATMOS_FORMATS = ["Marcs", "Stagger", "Stagger2", "Co5bold", "M3DIS", "Text", "save"]
STORAGE = "storage"        # the pseudo grid id of the storage folder's picker row


def _new_id(existing):
    """A short stable id (``g1``, ``g2``, …) not already in ``existing``."""
    for n in itertools.count(1):
        cand = f"g{n}"
        if cand not in existing:
            return cand


class GridManager:
    """Owns the grid page: storage settings, catalogue grids, custom grids, the jobs.

    ``on_change`` (optional) is called after any change so the panel can refresh its grid
    dropdown / file-existence warnings. ``default_storage`` is where grids are stored
    until the user picks a folder; ``ctx`` supplies ``submit_io`` / ``schedule`` /
    ``notify`` (without it the page renders but cannot download).
    """

    def __init__(self, default_pattern, on_change=None, default_storage=None, ctx=None):
        self.default_pattern = default_pattern
        self.on_change = on_change
        self.ctx = ctx
        self.default_storage = Path(default_storage).expanduser() if default_storage \
            else Path("~/atmos").expanduser()
        prefs = store.get_section("synthesis")
        # records with a catalogue id were the auto-registered bundled MARCS grid of the
        # earlier page; the catalogue row covers that folder now
        self.custom = [dict(g) for g in prefs.get("atmos_grids", [])
                       if g.get("id") not in grid_catalog.CATALOG]
        if len(self.custom) != len(prefs.get("atmos_grids", [])):
            store.update_pref("synthesis", "atmos_grids", [dict(g) for g in self.custom])
        # the settings are held here (only this page edits them), so the registry other
        # tools read never touches the store
        self._storage = (prefs.get("grid_storage") or "").strip()
        self._snapshots = str(prefs.get("grid_snapshots") or DEFAULT_SNAPSHOTS)
        self._editing = None      # id of the grid whose folder the picker is choosing
        self._open = None         # catalogue id whose model list is open
        self._listings = {}       # catalogue id -> Listing, once listed this session
        self._failed = {}         # catalogue id -> why its last listing failed
        self._running = {}        # catalogue id -> the job's status text while it runs
        self._status = ""
        # the model table's order, shared by every grid: stellar parameters, name, held
        self._sort = SortState({"teff": True, "logg": True, "feh": True, "name": True,
                                "held": False}, "teff")

        self.list_bridge = html_bridge("grid-list", name="grid-list-render")
        self.status_bridge = html_bridge("grid-status")
        self.open_bridge = class_bridge("#grid-picker", "open")
        self.fields_input = TextInput(value="", visible=False, name="grid-fields-input")
        self.actions = TextInput(value="", visible=False, name="grid-actions")
        self.fields_input.on_change("value", self._on_field)
        self.actions.on_change("value", self._on_action)

        self.picker = FilePicker(self._on_folder_pick, container_id="grid-results",
                                 start_dir=str(self.storage), pick_dirs=True)
        self.picker.layout.name = "grid_picker"

    # --- doc wiring -------------------------------------------------------
    @property
    def models(self):
        return [self.list_bridge, self.status_bridge, self.open_bridge, self.fields_input,
                self.actions]

    @property
    def results_name(self):
        return self.picker.results_name

    # --- settings ---------------------------------------------------------
    @property
    def storage(self):
        return Path(self._storage).expanduser() if self._storage else self.default_storage

    @property
    def snapshots(self):
        """The snapshots per model every download fetches and every catalogue grid's
        record selects."""
        return self._snapshots

    # --- the registry other tools read -----------------------------------
    def _catalog_records(self):
        out = []
        for spec in grid_catalog.CATALOG.values():
            if grid_catalog.models_on_disk(spec, self.storage):
                out.append(grid_catalog.grid_record(
                    spec, grid_catalog.grid_folder(spec, self.storage), self.snapshots))
        return out

    @property
    def grids(self):
        """Every registered grid record: the catalogue grids on disk, then the custom ones."""
        return self._catalog_records() + [dict(g) for g in self.custom]

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
        """Names of custom grids whose configured folder does not exist (empty ignored)."""
        out = []
        for g in self.custom:
            folder = (g.get("folder") or "").strip()
            if folder and not Path(folder).expanduser().is_dir():
                out.append(g.get("name") or g["id"])
        return out

    # --- rendering --------------------------------------------------------
    def list_html(self):
        if self._open is not None:
            body = self._models_html(self._open)
        else:
            body = (self._settings_html() + self._catalog_html()
                    + '<div class="section-header">Custom grids</div>'
                    + "".join(self._card_html(g) for g in self.custom)
                    + '<button class="text-btn grid-add"><i class="ti ti-plus"></i> '
                    'Add custom grid</button>' + self._ref_html())
        return (f'<div id="grid-status" class="field-static grid-status">'
                f'{html.escape(self._status)}</div>' + body)

    def render(self):
        self.list_bridge.value = self.list_html()

    def _settings_html(self):
        storage = str(self.storage)
        warn = " field-warn-input" if not self.storage.is_dir() else ""
        return (
            '<div class="section-header">Storage</div>'
            '<div class="field-group"><label>Storage folder</label>'
            f'<div class="path-row grid-folder-row" data-grid-id="{STORAGE}">'
            f'<input type="text" class="field-input grid-folder{warn}" '
            f'data-grid-id="{STORAGE}" data-grid-key="folder" value="{html.escape(storage)}" '
            'placeholder="Folder holding every downloaded grid" title="Click to browse">'
            '</div></div>'
            '<div class="field-group grid-count"><label>Snapshots per model</label>'
            f'<input type="text" class="field-input" data-grid-id="{STORAGE}" '
            f'data-grid-key="snapshots" value="{html.escape(self.snapshots)}" '
            'title="How many snapshots of each 3D model downloads fetch and runs use; '
            'raising it later adds snapshots to what is on disk"></div>')

    def _catalog_html(self):
        rows = ['<div class="section-header">Grids</div>']
        for spec in grid_catalog.CATALOG.values():
            gid = html.escape(spec.id)
            on_disk = grid_catalog.models_on_disk(spec, self.storage)
            listing = self._listings.get(spec.id)
            if spec.id in self._running:
                state = self._running[spec.id]
            elif not on_disk:
                state = "not downloaded"
            elif listing is not None:
                state = f"{len(on_disk)} of {len(listing.models)} models on disk"
            else:
                state = f"{len(on_disk)} models on disk"
            busy = " disabled" if spec.id in self._running else ""
            edit = (f'<button class="icon-btn grid-edit" data-grid-id="{gid}" '
                    f'title="Models"><i class="ti ti-list"></i></button>' if spec.snapshots
                    else "")
            delete = (f'<button class="icon-btn grid-del-catalog" data-grid-id="{gid}" '
                      f'data-name="{html.escape(spec.name)}" title="Delete the downloaded '
                      f'files"{busy}><i class="ti ti-trash"></i></button>' if on_disk else "")
            rows.append(
                f'<div class="grid-row" data-grid-id="{gid}">'
                f'<div class="grid-row-meta"><span class="grid-row-name" '
                f'title="{html.escape(spec.citation)}">{html.escape(spec.name)}</span>'
                f'<span class="grid-row-state">{html.escape(spec.format)} · '
                f'{html.escape(state)}</span></div>'
                f'<button class="icon-btn grid-fetch-all" data-grid-id="{gid}" '
                f'title="Download all models"{busy}><i class="ti ti-download"></i></button>'
                f'{edit}{delete}</div>')
        return "".join(rows)

    def _models_html(self, gid):
        spec = grid_catalog.CATALOG[gid]
        listing = self._listings.get(gid)
        head = (f'<div class="grid-card-head"><button class="icon-btn grid-back" '
                f'title="Back to the grids"><i class="ti ti-arrow-left"></i></button>'
                f'<span class="field-static grid-row-name">{html.escape(spec.name)}</span>'
                f'<a class="hint" href="{html.escape(spec.page)}" target="_blank" '
                f'rel="noopener">{html.escape(spec.page)}</a></div>')
        if listing is None:
            if gid in self._failed:
                return (head + f'<div class="field-static">Could not list the grid: '
                        f'{html.escape(self._failed[gid])}</div>'
                        f'<button class="text-btn grid-edit" data-grid-id="{gid}">'
                        '<i class="ti ti-refresh"></i> Retry</button>')
            return head + '<div class="field-static">Listing the grid…</div>'
        folder = grid_catalog.grid_folder(spec, self.storage)
        busy = " disabled" if gid in self._running else ""
        held = {m.name: len(m.on_disk(folder)) for m in listing.models}
        keys = {"teff": lambda m: m.teff, "logg": lambda m: m.logg, "feh": lambda m: m.feh,
                "name": lambda m: m.name, "held": lambda m: held[m.name]}
        order = self._sort.order(listing.models, keys[self._sort.col])
        rows = ['<div class="grid-models"><table class="linelist-table"><thead><tr>'
                + sort_headers(self._sort, [("Model", "name"), ("T<sub>eff</sub>", "teff"),
                                            ("log g", "logg"), ("[Fe/H]", "feh"),
                                            ("Snapshots", "held")])
                + '<th></th></tr></thead><tbody>']
        for m in (listing.models[i] for i in order):
            name = html.escape(m.name)
            fetch = (f'<button class="icon-btn grid-fetch-model" data-grid-id="{gid}" '
                     f'data-model="{name}" title="Download"{busy}>'
                     f'<i class="ti ti-download"></i></button>')
            delete = (f'<button class="icon-btn grid-del-model" data-grid-id="{gid}" '
                      f'data-model="{name}" data-name="{name}" title="Delete"{busy}>'
                      f'<i class="ti ti-trash"></i></button>' if held[m.name] else "")
            rows.append(
                f'<tr><td class="line-ion" title="{name}">{name}</td><td>{m.teff:.0f}</td>'
                f'<td>{m.logg:.2f}</td><td>{m.feh:+.2f}</td>'
                f'<td title="on disk / available">{held[m.name]} / {len(m.available)}</td>'
                f'<td>{fetch}{delete}</td></tr>')
        rows.append('</tbody></table></div>')
        return head + "".join(rows)

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
            # the folder row: focusing the field swaps in the shared folder picker (the
            # page moves #grid-picker here), whose listing drops down under it
            f'<div class="path-row grid-folder-row" data-grid-id="{gid}">'
            f'<input type="text" class="{folder_cls}" data-grid-id="{gid}" '
            f'data-grid-key="folder" value="{html.escape(folder)}" '
            f'placeholder="Atmosphere folder" title="Click to browse">'
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
                f'title="How many snapshots (in Halton order), \'all\', or explicit '
                f'snapshot numbers represent each node; node spectra are their mean">'
                f'<input type="text" class="field-input grid-resolution" data-grid-id="{gid}" '
                f'data-grid-key="resolution" value="{html.escape(resolution)}" '
                f'placeholder="Mesh nx,ny,nz (blank: Engine preferences)">')
        return "".join(out)

    def _ref_html(self):
        """The '1D LTE reference grid' selector, shown when grids exist.

        Names the registered grid used as the LTE baseline by the NLTE-correction modal.
        Hidden entirely when no grids are registered (there is nothing to choose)."""
        names = self.names
        if not names:
            return ""
        ref = self.ref_grid_id
        opts = "".join(
            f'<option value="{html.escape(gid)}"{" selected" if gid == ref else ""}>'
            f'{html.escape(name)}</option>' for gid, name in names)
        return ('<div class="section-header">Reference</div>'
                '<div class="field-group"><label>1D LTE reference grid</label>'
                f'<select class="field-input grid-ref-select">{opts}</select></div>')

    # --- state ------------------------------------------------------------
    def _persist(self):
        store.update_pref("synthesis", "atmos_grids", [dict(g) for g in self.custom])
        self._changed()

    def _changed(self):
        if self.on_change:
            self.on_change()

    def _set_status(self, text):
        self._status = text
        self.status_bridge.value = html.escape(text)

    def _custom(self, gid):
        return next((g for g in self.custom if g.get("id") == gid), None)

    def _on_field(self, attr, old, new):
        if not new:
            return
        msg = new.rsplit(" ", 1)[0]              # strip the trailing nonce
        ident, _, value = msg.partition("=")     # value may contain spaces / "{...}"
        gid, _, key = ident.partition("|")
        if gid == STORAGE:
            self._on_setting(key, value)
            return
        g = self._custom(gid)
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

    def _on_setting(self, key, value):
        value = value.strip()
        if key == "folder":
            self._storage = value
            store.update_pref("synthesis", "grid_storage", value)
            self._changed()
            self.render()
        elif key == "snapshots":
            if not value.isdigit() or int(value) < 1:
                self.render()       # restore the stored count
                return
            self._snapshots = value
            store.update_pref("synthesis", "grid_snapshots", int(value))
            self._changed()

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]             # strip the trailing nonce
        gid, _, model = arg.partition(" ")
        catalog = gid in grid_catalog.CATALOG
        if kind == "add":
            gid = _new_id({g["id"] for g in self.grids})
            self.custom.append({"id": gid, "name": "New grid", "folder": "",
                                "pattern": self.default_pattern, "format": ATMOS_FORMATS[0]})
            self._persist()
            self.render()
        elif kind == "del":
            self.custom = [g for g in self.custom if g.get("id") != arg]
            self._persist()
            self.render()  # ref_grid_id self-heals to a live grid via _ref_html on render
        elif kind == "setref":
            if self._find(arg) is not None:
                store.update_pref("synthesis", "ref_grid_id", arg)
                self._changed()
        elif kind == "browse":
            self._editing = arg
            if arg == STORAGE:
                start = str(self.storage)
            else:
                g = self._custom(arg)
                start = (g.get("folder") or "").strip() if g else ""
            # the field shows the folder being browsed, its contents drop down under it
            self.picker.preset(f"{Path(start).expanduser()}/" if start else "")
            pulse(self.open_bridge, "on")
        elif kind == "edit" and catalog:
            self._open = gid
            self._failed.pop(gid, None)
            self.render()
            if gid not in self._listings:
                self._list(gid)
        elif kind == "back":
            self._open = None
            self.render()
        elif kind == "sort":
            if self._sort.toggle(arg):
                self.render()
        elif kind == "fetch-all" and catalog:
            self._fetch(gid, None)
        elif kind == "fetch-model" and catalog and model:
            self._fetch(gid, model)
        elif kind == "del-catalog" and catalog and gid not in self._running:
            grid_catalog.delete_grid(grid_catalog.CATALOG[gid], self.storage)
            self._set_status(f"Deleted {grid_catalog.CATALOG[gid].name}")
            self._changed()
            self.render()
        elif kind == "del-model" and catalog and model and gid not in self._running:
            grid_catalog.delete_model(grid_catalog.CATALOG[gid], self.storage, model)
            self._set_status(f"Deleted {model}")
            self._changed()
            self.render()

    def _on_folder_pick(self, filepath):
        # the listing's "use this folder" row selects a directory; a file picked instead
        # stands for its parent directory
        p = Path(filepath)
        self._adopt_folder(p if p.is_dir() else p.parent)

    def _adopt_folder(self, folder):
        if folder is None:
            return
        if self._editing == STORAGE:
            self._on_setting("folder", str(folder))
        else:
            g = self._custom(self._editing)
            if g is None:
                return
            g["folder"] = str(folder)
            self._persist()
            self.render()
        pulse(self.open_bridge, "off")

    # --- jobs -------------------------------------------------------------
    def _list(self, gid, then=None):
        """List catalogue grid ``gid`` on the IO pool, then call ``then(listing)``."""
        if self.ctx is None:
            return
        spec = grid_catalog.CATALOG[gid]

        def work(token):
            try:
                return grid_catalog.list_models(gid, cancel=token)
            except CancelledError:
                self.ctx.schedule(lambda: self._list_failed(gid, "cancelled"))
                raise
            except Exception as exc:
                error = exc
                self.ctx.schedule(lambda: self._list_failed(gid, error))
                raise

        def done(listing):
            self._listings[gid] = listing
            self.render()
            if then is not None:
                then(listing)
        self.ctx.submit_io(f"list {spec.name}", work, done)

    def _list_failed(self, gid, exc):
        """The listing job raised: the page says why instead of waiting; the Processes
        row keeps the failure too."""
        self._failed[gid] = str(exc)
        self._set_status(f"Could not list {grid_catalog.CATALOG[gid].name}: {exc}")
        self.render()

    def _fetch(self, gid, model):
        """Fetch every model of ``gid`` (``model`` None) or one, with the shared count."""
        if self.ctx is None or gid in self._running:
            return
        if gid not in self._listings:
            self._list(gid, lambda listing: self._fetch(gid, model))
            return
        spec = grid_catalog.CATALOG[gid]
        plan = grid_catalog.plan(self._listings[gid], grid_catalog.grid_folder(spec, self.storage),
                                 int(self.snapshots), None if model is None else [model])
        what = spec.name if model is None else model
        if not plan.pending:
            self._set_status(f"{what} is already on disk")
            return
        self._running[gid] = f"downloading {grid_catalog.human(plan.pending_bytes)}"
        self._set_status(f"Downloading {what}: {grid_catalog.human(plan.pending_bytes)}")
        self.render()
        last = [0.0]

        def progress(done, total, label):
            now = time.monotonic()
            if now - last[0] < 0.5:
                return
            last[0] = now
            pct = 100 * done / total if total else 0
            text = (f"Downloading {what}: {pct:.0f}% · {grid_catalog.human(done)} of "
                    f"{grid_catalog.human(total)} · {label}")
            self.ctx.schedule(lambda: self._set_status(text))

        def work(token):
            try:
                return grid_catalog.fetch(plan, progress=progress, cancel=token)
            finally:
                self.ctx.schedule(lambda: self._settled(gid))

        def done(record):
            self._set_status(f"Downloaded {what}")
            self.ctx.notify(f"{what} downloaded")
        if self.ctx.submit_io(f"download {what}", work, done) is None:
            self._settled(gid)

    def _settled(self, gid):
        """The download of ``gid`` ended (finished, failed or cancelled): the rows show
        what is on disk now and the registry may have gained the grid."""
        self._running.pop(gid, None)
        self._changed()
        self.render()
