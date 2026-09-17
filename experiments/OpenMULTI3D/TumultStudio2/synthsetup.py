"""A reusable *synthesis setup* component: the atmosphere + stellar-parameter + RT-mode
cluster shared by any tool that builds a :class:`~tumult.synthesis.SynthRequest`.

It is the cluster of the synthesis modal (:class:`synthmodal.SynthModal`): grid selector +
Interpolate/Pick toggle + Format + atmosphere picker with filename→param auto-fill, plus
the Preset dropdown, the Teff/logg/[Fe/H]/vmic, abundances and LTE/NLTE forms. The
picker (:class:`atmospicker.AtmosPicker`) lists a registered grid's parsed models and
browses a folder for the Individual atmosphere. The Preset dropdown lists the reference
stars (:mod:`tumult.refstars`): it fills the four parameter fields from an entry and shows
*Custom* whenever the fields hold anything else, so the two stay in step in both
directions.

Parameterised by a ``prefix`` (DOM-id / form-name / action-channel namespace) and a store
``section`` (where this instance's per-target state lives), so several independent setups
can coexist — e.g. the synthesis modal and the Abundance tool's own modal.

The app-global **grid registry** (folders/patterns/formats) is *not* duplicated: a shared
:class:`gridmanager.GridManager` is passed in and used read-only; only the *selection*
(``atmos_grid_id``, ``atmos_mode``, ``atmos_file``, …) is per-section. NLTE model atoms are
likewise treated as a global, physics-level choice: instead of spawning a second
``AtomPicker``, an ``atoms_provider`` callable supplies the current selection for
:meth:`context` and a read-only display.

A host running a line-mask analysis sets :attr:`SynthSetup.element`; the mode form then
offers **Include background lines** (:attr:`SynthSetup.background_lines`), whose off state
the host passes on as ``background_lines=False`` to restrict the syntheses to that element.
In NLTE the synthesis modal offers **Use line list** (:attr:`SynthSetup.use_linelist`), off
meaning the model atoms' transitions alone.

``context()`` returns the :class:`~tumult.synthesis.SynthRequest` (without windows; the
caller adds its own) that both the synthesis and the fit mode of the modal start from, or
``None`` with a ``reason`` when the setup does not resolve to a supported request. NLTE
with model atoms is requested departures (``line_grid=False``: the atoms' own quadratures),
and saved departures are loaded as :class:`~tumult.departures.Departures`.
"""
from __future__ import annotations

import html
from pathlib import Path

from bokeh.models import TextInput

from tumult import refstars
from tumult.synthesis import (Abundances, Atmosphere, AtmosphereGrid, DepartureRequest,
                              Departures, Snapshot, Star, SynthRequest)
from tumult.synthesis.snapshot_grid import is_snapshot_grid

import store
from tumult.synthesis import workflows as synth
from atmospicker import AtmosPicker
from controls import Field, Form
from gridmanager import ATMOS_FORMATS
from htmlbridge import class_bridge, html_bridge, parse_action, pulse
from synthprefs import (INDIVIDUAL_GRID, REPO_DIR, RUN_DEFAULTS, atmosphere_grid, load_prefs,
                        missing_prefs, synth_config)

# keys copied when a setup inherits another section's configuration (the four stellar axes
# plus the abundance/mode/atmosphere selection); the grid id is handled separately because the
# synthesis modal stores it under "synthesis", not its run section.
_INHERIT_KEYS = ("teff", "logg", "feh", "vmic", "abundances", "relative_abundances", "mode",
                 "departures", "atmos_mode", "atmos_file",
                 "atmos_individual_format")

_SNAPSHOT_FORMATS = {"stagger", "stagger2", "co5bold", "m3dis", "must", "multi"}


def model_is_3d(atmos_format, path=None, root=None):
    """Infer dimensionality from the selected format, inspecting saved meshes as needed."""
    fmt = str(atmos_format or "").strip().lower()
    if fmt in _SNAPSHOT_FORMATS:
        return True
    if fmt != "save" or not path:
        return False
    model = Path(path).expanduser()
    if not model.is_absolute() and root is not None:
        model = Path(root).expanduser() / model
    try:
        nx, ny, _ = (int(value) for value in
                     (model / "atmos_mesh.txt").read_text().splitlines()[0].split()[:3])
        return nx * ny > 1
    except (OSError, ValueError, IndexError):
        return False


class SynthSetup:
    """Atmosphere + stellar-parameter + RT-mode configuration for one consumer.

    Parameters
    ----------
    prefix : str
        Namespace for DOM ids / form names / the grid-action channel (e.g. ``"abund"`` →
        ``#abund-grid``, ``abund-grid-actions``, forms ``abundparams``/``abundbase``/
        ``abundmode``). The atmosphere selects carry ``data-atmos="{prefix}"``, which is
        how :data:`ATMOS_SCRIPT` finds the channel; choose CSS-class-safe values.
    section : str
        Store section holding this setup's per-target state (e.g. ``"abund_run"``).
    grids : gridmanager.GridManager
        The shared, app-global grid registry (used read-only).
    atoms_provider : callable | None
        ``() -> list[str]`` returning the current NLTE model atoms (shared, e.g. the
        synthesis modal's selection). Used by :meth:`context` and the read-only display.
    notify : callable | None
        Footer message sink.
    on_mode_change : callable | None
        ``(is_nlte: bool) -> None`` fired when the RT mode toggles (so the host can refresh
        the read-only atoms display).
    on_base_change : callable | None
        ``(key, value) -> None`` fired after an abundances / relative-abundances edit is
        persisted (so a host can sync dependent displays, e.g. the Fit tab's rows).
    on_run_change : callable | None
        ``(key, value) -> None`` fired after a stellar-parameter edit is persisted (so a host
        that states the setup elsewhere, e.g. the Abundance object panel, can restate it).
    grid_section : str | None
        Store section holding the grid *selection* (``atmos_grid_id``). Defaults to
        ``section``; the synthesis modal historically keeps it under ``"synthesis"``.
    form_names : dict | None
        Overrides for the three form names (keys ``"run"``/``"base"``/``"mode"``), for
        hosts whose page slots predate the ``{prefix}params``… convention.
    """

    def __init__(self, prefix, section, grids, *, atoms_provider=None, notify=None,
                 on_mode_change=None, on_base_change=None, on_atmos_change=None,
                 on_run_change=None, grid_section=None, form_names=None):
        self.allow_departures = prefix == "synth"
        self.allow_atom_only = prefix == "synth"
        # the element a line-mask analysis diagnoses; while set, the mode form offers to
        # synthesise that element's transitions alone (``only_elements``)
        self.element = ""
        self.prefix = prefix
        self.section = section
        self.grid_section = grid_section or section
        self._sections = {name: store.get_section(name) for name in {section, self.grid_section}}
        self.grids = grids
        self.atoms_provider = atoms_provider or (lambda: [])
        self.notify = notify or (lambda *a, **k: None)
        self.on_mode_change = on_mode_change
        self.on_base_change = on_base_change
        self.on_atmos_change = on_atmos_change    # fired when the atmosphere source changes
        self.on_run_change = on_run_change        # fired when a stellar parameter is edited
        self.reason = None                 # why the last context() returned None

        names = {"run": f"{prefix}params", "base": f"{prefix}base",
                 "mode": f"{prefix}mode", **(form_names or {})}
        vals = {**RUN_DEFAULTS, **store.get_section(section)}
        self.run_form = Form(names["run"], [
            # the reference-star catalogue as presets, full width like the Grid select
            # above it since the entry names are long
            Field("refstar", "Preset", kind="select", value=refstars.CUSTOM_KEY,
                  options=self._refstar_options()),
            *(Field(key, label, kind="number", value=vals[key], step=step, inline=True,
                    min=0 if key == "vmic" else None)
              for key, label, step, _ in synth.STELLAR_PARAMS),
        ], on_change=self._save_run)
        self._sync_refstar(render=False)
        self.base_form = Form(names["base"], [
            Field("abundances", "Abundance overrides", kind="text", value=vals["abundances"],
                  placeholder="e.g. Ba=2.3, Eu=0.5"),
            Field("relative_abundances", "Relative abundances", kind="checkbox",
                  value=vals["relative_abundances"]),
        ], on_change=self._save_base)
        self.mode_form = Form(names["mode"], [
            Field("use_linelist", "Use line list", kind="checkbox", value=vals.get("use_linelist", True)),
            Field("background_lines", "Include background lines", kind="checkbox",
                  value=vals.get("background_lines", True)),
            Field("mode", "Radiative transfer", kind="select", value=vals["mode"],
                  options=["LTE", "NLTE"], inline=True),
            Field("departures", "Saved departures (optional)", kind="text",
                  value=vals.get("departures", ""),
                  placeholder="path to departure run"),
        ], on_change=self._save_mode)

        # atmosphere selection UI: a server-rendered block (grid select + source/format + file
        # chip) pushed via an html bridge, edited through a client→server action channel.
        self.grid_bridge = html_bridge(f"{prefix}-grid", name=f"{prefix}-grid-render")
        self.grid_actions = TextInput(value="", visible=False, name=f"{prefix}-grid-actions")
        self.grid_actions.on_change("value", self._on_grid_action)
        # `open` reveals the always-on path field when a concrete file is in play; the
        # dropdown result list is toggled client-side via `results-open` (on focus), and
        # collapsed from the server (`file_results_bridge`) once a file is picked.
        self.file_open_bridge = class_bridge(f"#{prefix}-file-picker", "open")
        self.file_results_bridge = class_bridge(f"#{prefix}-file-picker", "results-open")
        self.atmos_picker = AtmosPicker(self._on_pick_atmos,
                                        container_id=f"{prefix}-file-results",
                                        root=self._individual_dir())
        self.atmos_picker.layout.name = f"{prefix}_file_picker"
        # list the selected grid's models (or the Individual folder) and seed the single
        # field with the persisted selection so it shows the current file
        self._nodes_stale = False
        self._reset_atmos_picker()
        self._seed_atmos_query()

    # --- models / rendering -----------------------------------------------------------
    @property
    def models(self):
        return [*self.run_form.models, *self.base_form.models, *self.mode_form.models,
                self.grid_bridge, self.grid_actions, self.file_open_bridge,
                self.file_results_bridge]

    def render(self):
        """Push every server-rendered surface (forms + grid) to the client."""
        self.run_form.render()
        self.base_form.render()
        self._sync_mode_fields()
        self.render_grid_ui()

    def _sync_mode_fields(self):
        """Show the mode form's conditional fields for the current mode and job, and repaint."""
        fields = self.mode_form.fields
        fields["use_linelist"].hidden = not self.allow_atom_only or not self.is_nlte
        fields["background_lines"].hidden = not self.element
        fields["departures"].hidden = not self.is_nlte or not self.allow_departures
        self.mode_form.render()

    def render_grid_ui(self):
        self.grid_bridge.value = self.grid_ui_html()
        # reveal the always-on path field exactly when a concrete file is in play
        pulse(self.file_open_bridge, "on" if self.show_file() else "off")
        # every atmosphere-source change (grid/mode/file) routes through here; let listeners
        # (e.g. the Fit tab) react to the new source.
        if self.on_atmos_change:
            self.on_atmos_change()

    def show_file(self):
        """True when a concrete atmosphere *file* is in play — Individual, or a registered
        grid in 'Pick model' mode — i.e. when the inline picker field should be shown."""
        if self.is_snapshot or self.is_individual:
            return True
        return self._sections[self.section].get("atmos_mode", "interp") == "pick"

    def _seed_atmos_query(self):
        """Show the currently-chosen file in the path field (empty when none / not picking):
        its model label on a grid, its path for an Individual atmosphere."""
        picked = (self._sections[self.section].get("atmos_file") or "").strip()
        if not self.show_file():
            picked = ""
        elif self._nodes_stale:
            self._reset_atmos_picker()     # the grid is listed once it is to be picked from
        self.atmos_picker.preset(self.atmos_picker.label_for(picked) or picked)

    # --- reference stars ----------------------------------------------------------
    @staticmethod
    def _refstar_options():
        return [(refstars.CUSTOM_KEY, "Custom")] + [(s.key, s.label) for s in refstars.stars()]

    def stellar_values(self):
        """``{teff, logg, feh, vmic}`` as the run form holds them."""
        return {k: self.run_form.value(k) for k, _, _, _ in synth.STELLAR_PARAMS}

    def _sync_refstar(self, render=True):
        """Point the dropdown at the entry whose parameters the fields hold, else Custom."""
        try:
            star = refstars.matching(**self.stellar_values())
        except (TypeError, ValueError):
            star = None
        key = star.key if star else refstars.CUSTOM_KEY
        if key != self.run_form.value("refstar"):
            self.run_form.fields["refstar"].value = key
            if render:
                self.run_form.render()

    def refresh_refstars(self):
        """The catalogue changed (Preferences → Reference stars): rebuild the presets."""
        self.run_form.fields["refstar"].options = self._refstar_options()
        self._sync_refstar(render=False)
        self.run_form.render()

    def set_stellar(self, values):
        """Write ``{teff, logg, feh, vmic}`` (any subset) into the form and this section,
        keeping the reference-star dropdown in step."""
        vals = {k: v for k, v in values.items() if k in self.stellar_values()}
        self.run_form.set_values(**vals)
        for k, v in vals.items():
            self._update(self.section, k, v)
        self._sync_refstar()

    # --- grid selection (per-section) -------------------------------------------------
    def selected_grid_id(self):
        """The selected grid id, validated against the live registry plus the Individual
        pseudo-grid (falls back to the first registered grid, else Individual)."""
        ids = [gid for gid, _ in self.grids.names] + [INDIVIDUAL_GRID]
        cur = self._sections[self.grid_section].get("atmos_grid_id")
        if cur not in ids:
            cur = ids[0]
            self._update(self.grid_section, "atmos_grid_id", cur)
        return cur

    @property
    def is_individual(self):
        return self.selected_grid_id() == INDIVIDUAL_GRID

    @property
    def is_nlte(self):
        return str(self.mode_form.value("mode")).upper() == "NLTE"

    @property
    def is_3d_format(self):
        """Whether the selected grid or individual format represents 3D atmospheres."""
        fmt = str((self.active_grid() or {}).get("format") or "").strip().lower()
        selected = (self._sections[self.section].get("atmos_file") or "").strip()
        root = str(load_prefs().get("dispatch_cwd") or REPO_DIR)
        return model_is_3d(fmt, selected, root)

    @property
    def is_snapshot(self):
        """Whether one 3D snapshot file is synthesised directly. A registered snapshot
        grid (see :mod:`tumult.synthesis.snapshot_grid`) may instead interpolate between
        its nodes like a 1D grid, so 'Interpolate' on such a grid is not a snapshot run."""
        if not self.is_3d_format:
            return False
        if self.is_individual:
            return True
        grid = self.selected_grid() or {}
        return (not is_snapshot_grid(grid)
                or self._sections[self.section].get("atmos_mode", "interp") == "pick")

    def selected_grid(self):
        """The selected *registered* grid dict, or ``None`` (no grids, or Individual)."""
        return None if self.is_individual else self.grids._find(self.selected_grid_id())

    def active_grid(self):
        """The grid dict to synthesise with: a registered grid, or a synthetic dict for the
        Individual pseudo-grid (no folder/pattern; format from the Format selector)."""
        if self.is_individual:
            fmt = self._sections[self.section].get("atmos_individual_format") or ATMOS_FORMATS[0]
            return {"name": "Individual", "folder": "", "pattern": "", "format": fmt}
        return self.grids._find(self.selected_grid_id())

    def atmos_file(self):
        """The concrete atmosphere file path when one is in play (Individual or a picked grid
        member), else ``None`` — the caller then interpolates the grid."""
        if not self.show_file():
            return None
        return (self._sections[self.section].get("atmos_file") or "").strip() or None

    def _individual_dir(self):
        cwd = Path(str(load_prefs().get("dispatch_cwd") or REPO_DIR)).expanduser()
        d = cwd / "input" / "atmos"
        return d if d.is_dir() else cwd

    def _reset_atmos_picker(self):
        """Browse the Individual folder, or list the selected grid's parsed models. A grid
        is only scanned while its models are to be picked from (:meth:`show_file`); until
        then the list is left empty and marked stale, so an interpolated setup never pays
        for the scan."""
        self._nodes_stale = False
        if self.is_individual:
            self.atmos_picker.set_root(self._individual_dir())
            return
        if not self.show_file():
            self.atmos_picker.set_nodes([])
            self._nodes_stale = True
            return
        grid = self.selected_grid()
        if grid is None:
            self.atmos_picker.set_nodes([], "No atmosphere grid registered")
            return
        root = str(load_prefs().get("dispatch_cwd") or REPO_DIR)
        self.atmos_picker.set_grid(grid, root)

    # --- atmosphere UI rendering ------------------------------------------------------
    def grid_ui_html(self):
        """Grid selector (registered grids + 'Individual atmosphere'); then a Source toggle
        (Interpolate / Pick model) for a registered grid, or a Format selector for Individual.
        The chosen model itself lives in the always-on picker field below this block (revealed
        by :meth:`show_file`), so no chip/browse button is rendered here. Each select carries
        ``data-atmos`` (this setup's prefix) and ``data-atmos-key`` (``grid``/``mode``/
        ``format``), the attributes :data:`ATMOS_SCRIPT` routes a change by."""
        p = self.prefix
        sel = self.selected_grid_id()
        run = self._sections[self.section]
        mode = run.get("atmos_mode", "interp")
        grid_opts = "".join(
            f'<option value="{html.escape(gid)}"{" selected" if gid == sel else ""}>'
            f'{html.escape(name)}</option>' for gid, name in self.grids.names)
        grid_opts += (f'<option value="{INDIVIDUAL_GRID}"'
                      f'{" selected" if sel == INDIVIDUAL_GRID else ""}>'
                      'Individual atmosphere</option>')
        out = ['<div class="field-group"><label>Grid</label>'
               f'<select class="field-input" data-atmos="{p}" data-atmos-key="grid">'
               f'{grid_opts}</select></div>']

        if sel == INDIVIDUAL_GRID:
            fmt = run.get("atmos_individual_format") or ATMOS_FORMATS[0]
            fmt_opts = "".join(
                f'<option value="{html.escape(f)}"{" selected" if f == fmt else ""}>'
                f'{html.escape(f)}</option>' for f in ATMOS_FORMATS)
            out.append('<div class="field-group field-inline"><label>Format</label>'
                       f'<select class="field-input" data-atmos="{p}" data-atmos-key="format">'
                       f'{fmt_opts}</select></div>')
        else:
            selected = self.selected_grid() or {}
            if model_is_3d(selected.get("format")) and not is_snapshot_grid(selected):
                out.append('<div class="field-group field-inline"><label>Source</label>'
                           '<div class="field-static">Snapshot file</div></div>')
            else:
                mode_opts = "".join(
                    f'<option value="{v}"{" selected" if v == mode else ""}>{lbl}</option>'
                    for v, lbl in (("interp", "Interpolate"), ("pick", "Pick model")))
                out.append('<div class="field-group field-inline"><label>Source</label>'
                           f'<select class="field-input" data-atmos="{p}" '
                           f'data-atmos-key="mode">{mode_opts}</select></div>')
        return "".join(out)

    # --- action handling --------------------------------------------------------------
    def _on_grid_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "grid":
            self._update(self.grid_section, "atmos_grid_id", arg)
            self._update(self.section, "atmos_file", "")  # drop a stale cross-grid file
            self._reset_atmos_picker()
            self._seed_atmos_query()
            self.render_grid_ui()
        elif kind == "mode":
            self._update(self.section, "atmos_mode", arg)
            self._seed_atmos_query()   # entering 'pick' shows the stored file (if any)
            self.render_grid_ui()
        elif kind == "format":
            self._update(self.section, "atmos_individual_format", arg)
        elif kind == "clearfile":
            self._update(self.section, "atmos_file", "")
            self._seed_atmos_query()
            self.render_grid_ui()

    def _on_pick_atmos(self, filepath):
        """A direct atmosphere file was picked: remember it and auto-fill the param fields
        from the filename via the grid's pattern (editable; non-matching files leave the
        fields as-is)."""
        path = str(filepath)
        self._update(self.section, "atmos_file", path)
        # the single field now shows the chosen model (its label on a grid, else its path)
        self.atmos_picker.query.value = self.atmos_picker.label_for(path) or path
        grid = self.selected_grid()
        if grid:
            # the pattern (the MARCS one when the record leaves it blank, as discovery
            # does) binds the path relative to the grid folder (a directory-per-node grid
            # names its parameters in the folder), or the bare name outside it
            name = Path(path).name
            folder = (grid.get("folder") or "").strip()
            if folder:
                try:
                    name = Path(path).expanduser().resolve().relative_to(
                        Path(folder).expanduser().resolve()).as_posix()
                except ValueError:
                    pass
            rec = synth.recover_params(grid.get("pattern") or synth.DEFAULT_PATTERN, name)
            if rec:
                teff, logg, feh, vmic = rec
                vals = {"teff": teff, "logg": logg, "feh": feh}
                if vmic is not None:
                    vals["vmic"] = vmic
                self.set_stellar(vals)
        self.render_grid_ui()
        pulse(self.file_results_bridge, "off")   # collapse the dropdown, keep the field

    # --- persistence ------------------------------------------------------------------
    def _save_run(self, key, value):
        if key == "refstar":
            # choosing an entry fills the four fields from it; Custom changes nothing
            star = refstars.by_key().get(str(value))
            if star is None:
                self._sync_refstar()
                return
            self.set_stellar(star.params)
            for k, v in star.params.items():
                if self.on_run_change:
                    self.on_run_change(k, v)
            return
        self._update(self.section, key, value)
        # typing over an entry's numbers is what Custom means
        self._sync_refstar()
        if self.on_run_change:
            self.on_run_change(key, value)

    def _save_base(self, key, value):
        self._update(self.section, key, value)
        if key == "relative_abundances":
            # convert fixed overrides so the physical abundance is preserved, anchored at
            # the synthesis [Fe/H].
            text = str(self.base_form.value("abundances") or "").strip()
            if text:
                feh = float(self.run_form.value("feh") or 0.0)
                solar = synth.abund_table(load_prefs())
                converted = synth.convert_abundances(text, to_relative=bool(value),
                                                     feh=feh, solar=solar)
                self.base_form.set_values(abundances=converted)
                self._update(self.section, "abundances", converted)
        if self.on_base_change:
            self.on_base_change(key, value)

    def _save_mode(self, key, value):
        self._update(self.section, key, value)
        self._sync_mode_fields()
        if self.on_mode_change:
            self.on_mode_change(self.is_nlte)

    # --- inherit / reload -------------------------------------------------------------
    def inherit_from(self, run_section="synth_run", grid_section="synthesis",
                     grid_key="atmos_grid_id", *, force=False):
        """Seed this setup from another section's configuration. By default only fills keys
        this section does not already have (first-open inheritance); ``force=True`` overwrites
        (the explicit 'Copy from the synthesis setup' action)."""
        dest = dict(self._sections[self.section])  # keeps the _seeded marker, if set
        src = store.get_section(run_section)
        for k in _INHERIT_KEYS:
            if k in src and (force or k not in dest):
                dest[k] = src[k]
        # always adopt the source's grid *selection* on inherit (the construction-time default
        # in this section must not suppress it); the synthesis modal stores the id under
        # ``grid_section``/``grid_key`` ("synthesis"/"atmos_grid_id"), not its run section.
        gid = store.get_section(grid_section).get(grid_key)
        if gid:
            dest["atmos_grid_id"] = gid
        self._sections[self.section] = dest
        self.reload()

    def _update(self, section, key, value):
        self._sections[section][key] = value

    def commit(self):
        """Remember the submitted setup without persisting unrelated preference keys."""
        self._sections[self.section].update(self.stellar_values())
        self._sections[self.section].update(self.base_form.values())
        self._sections[self.section].update(self.mode_form.values())
        with store.transaction() as data:
            prefs = data["preferences"]
            prefs.setdefault(self.section, {}).update(self._sections[self.section])
            if self.grid_section != self.section:
                prefs.setdefault(self.grid_section, {})["atmos_grid_id"] = self.selected_grid_id()

    def discard(self):
        """Restore the last submitted configuration after dismissing a draft."""
        self._sections = {name: store.get_section(name)
                          for name in {self.section, self.grid_section}}
        self.reload()

    def refresh_grids(self):
        """The grid registry changed: re-list the selected grid's models and repaint the
        selector (a grid that left the registry falls back in :meth:`selected_grid_id`)."""
        self._reset_atmos_picker()
        self._seed_atmos_query()
        self.render_grid_ui()

    def reload(self):
        """Re-seed the forms + atmosphere UI from the (possibly just-changed) store section."""
        vals = {**RUN_DEFAULTS, **self._sections[self.section]}
        self.run_form.set_values(**{k: vals[k] for k, _, _, _ in synth.STELLAR_PARAMS})
        self._sync_refstar()
        self.base_form.set_values(abundances=vals["abundances"],
                                  relative_abundances=vals["relative_abundances"])
        self.mode_form.set_values(mode=vals["mode"], use_linelist=vals.get("use_linelist", True),
                                  background_lines=vals.get("background_lines", True),
                                  departures=vals.get("departures", ""))
        self._reset_atmos_picker()
        self._seed_atmos_query()
        self.render_grid_ui()

    def load_request(self, request, use_linelist=True):
        """Show ``request``'s atmosphere, star, abundance overrides and RT mode in this
        setup, the inverse of :meth:`context`. It edits the unsaved state only, so a draft
        discard restores the last submitted setup. The configuration (paths, resolution)
        stays with Preferences. Returns why the atmosphere could not be selected, or
        ``None``."""
        sect = self._sections[self.section]
        sect.update({k: v for k, v in zip(("teff", "logg", "feh", "vmic"), request.star.tuple())
                     if v is not None})
        abundances = request.abundances
        sect["abundances"] = ", ".join(f"{el}={value:g}"
                                       for el, value in (abundances or {}).items())
        if abundances is not None:
            sect["relative_abundances"] = abundances.relative
        sect["mode"] = "NLTE" if request.is_nlte else "LTE"
        sect["departures"] = (str(request.departures.folder)
                              if isinstance(request.departures, Departures) else "")
        sect["use_linelist"] = bool(use_linelist)
        reason = self._load_atmosphere(request.atmosphere, request.config)
        self.reload()
        return reason

    def _load_atmosphere(self, atmosphere, config):
        """Select the registered grid ``atmosphere`` names (interpolated), the grid holding
        its file (picked), or the Individual atmosphere with its format."""
        sect = self._sections[self.section]
        grids = [(gid, self.grids._find(gid) or {}) for gid, _ in self.grids.names]
        if isinstance(atmosphere, AtmosphereGrid):
            folder = config.resolve(atmosphere.folder)
            gid = next((gid for gid, g in grids
                        if config.resolve(g.get("folder") or "") == folder), None)
            if gid is None:
                return (f"The atmosphere grid {atmosphere.folder} is not registered; the "
                        "grid selection was left as it is.")
            sect.update(atmos_mode="interp", atmos_file="")
        else:
            path = config.resolve(atmosphere.path)
            gid = next((gid for gid, g in grids if g.get("folder")
                        and path.is_relative_to(config.resolve(g["folder"]))),
                       INDIVIDUAL_GRID)
            sect.update(atmos_mode="pick", atmos_file=str(atmosphere.path))
            if gid == INDIVIDUAL_GRID:
                sect["atmos_individual_format"] = atmosphere.format
        self._sections[self.grid_section]["atmos_grid_id"] = gid
        return None

    # --- the payload ------------------------------------------------------------------
    @property
    def use_linelist(self):
        """False when an NLTE synthesis is to hold the model atoms' transitions alone."""
        return not (self.allow_atom_only and self.is_nlte
                    and not self.mode_form.value("use_linelist"))

    @property
    def background_lines(self):
        """False when a line-mask analysis synthesises its element's transitions alone."""
        return not (self.element and not self.mode_form.value("background_lines"))

    def context(self):
        """A :class:`~tumult.synthesis.SynthRequest` for the operations, or ``None`` (with
        :attr:`reason` set) when synthesis is not yet runnable or the setup does not form a
        valid request. The request carries no windows; the caller adds its own."""
        prefs = load_prefs()
        if missing_prefs(prefs):
            self.reason = "Synthesis paths are missing — check Preferences → Synthesis."
            return None
        grid = self.active_grid()
        if grid is None:
            self.reason = ("No atmosphere grid — download or register one in Preferences → Atmosphere grids, "
                           "or choose 'Individual atmosphere'.")
            return None
        stellar = self.stellar_values()
        atoms = list(self.atoms_provider())
        if self.is_nlte and not atoms:
            self.reason = "NLTE needs model atoms — add them in this calculation."
            return None
        sect = self._sections[self.section]
        atmos_file = None
        if self.is_snapshot or self.is_individual or sect.get("atmos_mode", "interp") == "pick":
            atmos_file = (sect.get("atmos_file") or "").strip() or None
            if atmos_file is None:
                self.reason = "Choose a grid to interpolate, or pick an atmosphere file."
                return None
        saved = str((self.mode_form.value("departures") if not self.mode_form.fields["departures"].hidden else "") or "").strip()
        try:
            departures = Departures.load(saved) if saved else None
        except (OSError, ValueError) as exc:
            self.reason = f"Cannot open saved departures: {exc}"
            return None
        if departures is None and self.is_nlte:
            departures = DepartureRequest(atoms, line_grid=False)
        reference = self.grids.reference_grid()
        preparation = atmosphere_grid(reference) if reference else None
        try:
            if self.is_snapshot:
                resolution = tuple(int(v) for v in (grid.get("resolution") or
                                                    [prefs.get(k) for k in ("nx", "ny", "nz")]))
                # active_grid() already carries the Individual selector's format or the
                # registered grid's; a snapshot picked from a Stagger grid must not fall back
                # to the Individual preference (which defaults to Marcs)
                atmosphere = Snapshot(atmos_file, grid.get("format"), resolution, preparation)
            elif atmos_file:
                atmosphere = Atmosphere(atmos_file, grid.get("format") or "Marcs", preparation)
            else:
                atmosphere = atmosphere_grid(grid, preparation)
            base = self.base_form.values()
            convention = "[X/Fe]" if base["relative_abundances"] else "A(X)"
            request = SynthRequest(
                config=synth_config(prefs), atmosphere=atmosphere,
                star=Star(*(float(stellar[k]) for k in ("teff", "logg", "feh", "vmic"))),
                abundances=Abundances.from_override(base["abundances"], convention),
                departures=departures)
        except (TypeError, ValueError) as exc:
            self.reason = str(exc)
            return None
        self.reason = None
        return request


# One delegated listener drives the atmosphere selects of every SynthSetup on the page
# (the synthesis modal, the abundance modal's Setup tab, the corrections modal). Emitted once
# through the ``atmos_script`` template variable, beside the Form script. The payload is the
# space-delimited ``"<key> <value> <nonce>"`` that :func:`htmlbridge.parse_action` reads, so
# a value keeps its spaces; the file picker's own path field is not part of this.
ATMOS_SCRIPT = r"""
(function () {
  var seq = 0;   // a repeated identical choice within one millisecond is still a new value
  document.addEventListener('change', function (ev) {
    var el = ev.target.closest('select[data-atmos]');
    if (!el) return;
    var m = Bokeh.documents[0].get_model_by_name(el.dataset.atmos + '-grid-actions');
    if (m) m.setv({ value: el.dataset.atmosKey + ' ' + el.value + ' ' + Date.now() + '.' + (++seq) });
  });
})();
"""
