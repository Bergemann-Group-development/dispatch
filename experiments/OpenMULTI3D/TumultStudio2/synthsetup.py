"""A reusable *synthesis setup* component: the atmosphere + stellar-parameter + RT-mode
cluster shared by any tool that needs to drive ``synth.compute``.

It is a generalisation of the cluster originally embedded in :class:`synthpanel.SynthPanel`
(grid selector + Interpolate/Pick toggle + Format + atmosphere file picker with filename→param
auto-fill, plus the Teff/logg/[Fe/H]/vmic, abundances and LTE/NLTE forms). Parameterised by a
``prefix`` (DOM-id / form-name / action-channel namespace) and a store ``section`` (where this
instance's per-target state lives), so several independent setups can coexist — e.g. the
Synthesize tab (``synthpanel``, untouched) and the Abundance tool's own Setup tab.

The app-global **grid registry** (folders/patterns/formats) is *not* duplicated: a shared
:class:`gridmanager.GridManager` is passed in and used read-only; only the *selection*
(``atmos_grid_id``, ``atmos_mode``, ``atmos_file``, …) is per-section. NLTE model atoms are
likewise treated as a global, physics-level choice: instead of spawning a second
``AtomPicker``, an ``atoms_provider`` callable supplies the current selection for
:meth:`context` and a read-only display.

``context()`` returns the ``{prefs, grid, run_vals, target, atmos_file}`` dict that
``synth.compute`` consumes (or ``None`` with a ``reason`` when no atmosphere resolves).
"""
from __future__ import annotations

import html
from pathlib import Path

from bokeh.models import TextInput

from tumult.synthesis import SynthRequest

import store
import synth
from controls import Field, Form
from filepicker import FilePicker
from gridmanager import ATMOS_FORMATS
from htmlbridge import class_bridge, html_bridge, parse_action, pulse
from synthprefs import INDIVIDUAL_GRID, REPO_DIR, RUN_DEFAULTS, load_prefs, missing_prefs

# keys copied when a setup inherits another section's configuration (the four stellar axes
# plus the abundance/mode/atmosphere selection); the grid id is handled separately because the
# Synthesize tab stores it under "synthesis", not its run section.
_INHERIT_KEYS = ("teff", "logg", "feh", "vmic", "abundances", "relative_abundances", "mode",
                 "atmos_mode", "atmos_file", "atmos_individual_format")


class SynthSetup:
    """Atmosphere + stellar-parameter + RT-mode configuration for one consumer.

    Parameters
    ----------
    prefix : str
        Namespace for DOM ids / form names / the grid-action channel (e.g. ``"abund"`` →
        ``#abund-grid``, ``abund-grid-actions``, forms ``abundparams``/``abundbase``/
        ``abundmode``). Choose CSS-class-safe values; the page JS keys on ``.{prefix}-…``.
    section : str
        Store section holding this setup's per-target state (e.g. ``"abund_run"``).
    grids : gridmanager.GridManager
        The shared, app-global grid registry (used read-only).
    atoms_provider : callable | None
        ``() -> list[str]`` returning the current NLTE model atoms (shared, e.g. the
        Synthesize tab's selection). Used by :meth:`context` and the read-only display.
    notify : callable | None
        Footer message sink.
    on_mode_change : callable | None
        ``(is_nlte: bool) -> None`` fired when the RT mode toggles (so the host can refresh
        the read-only atoms display).
    on_base_change : callable | None
        ``(key, value) -> None`` fired after an abundances / relative-abundances edit is
        persisted (so a host can sync dependent displays, e.g. the Fit tab's rows).
    grid_section : str | None
        Store section holding the grid *selection* (``atmos_grid_id``). Defaults to
        ``section``; the Synthesize tab historically keeps it under ``"synthesis"``.
    form_names : dict | None
        Overrides for the three form names (keys ``"run"``/``"base"``/``"mode"``), for
        hosts whose page slots predate the ``{prefix}params``… convention.
    """

    def __init__(self, prefix, section, grids, *, atoms_provider=None, notify=None,
                 on_mode_change=None, on_base_change=None, on_atmos_change=None,
                 grid_section=None, form_names=None):
        self.prefix = prefix
        self.section = section
        self.grid_section = grid_section or section
        self.grids = grids
        self.atoms_provider = atoms_provider or (lambda: [])
        self.notify = notify or (lambda *a, **k: None)
        self.on_mode_change = on_mode_change
        self.on_base_change = on_base_change
        self.on_atmos_change = on_atmos_change    # fired when the atmosphere source changes
        self.reason = None                 # why the last context() returned None

        names = {"run": f"{prefix}params", "base": f"{prefix}base",
                 "mode": f"{prefix}mode", **(form_names or {})}
        vals = {**RUN_DEFAULTS, **store.get_section(section)}
        self.run_form = Form(names["run"], [
            Field(key, label, kind="number", value=vals[key], step=step, inline=True,
                  min=0 if key == "vmic" else None)
            for key, label, step, _ in synth.STELLAR_PARAMS
        ], on_change=self._save_run)
        self.base_form = Form(names["base"], [
            Field("abundances", "Abundance overrides", kind="text", value=vals["abundances"],
                  placeholder="e.g. Ba=2.3, Eu=0.5"),
            Field("relative_abundances", "Relative abundances", kind="checkbox",
                  value=vals["relative_abundances"]),
        ], on_change=self._save_base)
        self.mode_form = Form(names["mode"], [
            Field("mode", "RT mode", kind="select", value=vals["mode"],
                  options=["LTE", "NLTE"], inline=True),
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
        self.atmos_picker = FilePicker(self._on_pick_atmos,
                                       container_id=f"{prefix}-file-results",
                                       root=self._atmos_picker_root())
        self.atmos_picker.layout.name = f"{prefix}_file_picker"
        # seed the single field with the persisted selection so it shows the current file
        self._seed_atmos_query()
        # read-only NLTE atoms note (shared with the Synthesize tab)
        self.atoms_ro_bridge = html_bridge(f"{prefix}-atoms-ro",
                                            name=f"{prefix}-atoms-ro-render")

    # --- models / rendering -----------------------------------------------------------
    @property
    def models(self):
        return [*self.run_form.models, *self.base_form.models, *self.mode_form.models,
                self.grid_bridge, self.grid_actions, self.file_open_bridge,
                self.file_results_bridge, self.atoms_ro_bridge]

    def render(self):
        """Push every server-rendered surface (forms + grid + atoms note) to the client."""
        self.run_form.render()
        self.base_form.render()
        self.mode_form.render()
        self.render_grid_ui()
        self.render_atoms_ro()

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
        grid in 'Pick file' mode — i.e. when the inline path field should be shown."""
        if self.is_individual:
            return True
        return store.get_section(self.section).get("atmos_mode", "interp") == "pick"

    def _seed_atmos_query(self):
        """Show the currently-chosen file in the path field (empty when none / not picking)."""
        picked = (store.get_section(self.section).get("atmos_file") or "").strip()
        self.atmos_picker.preset(picked if self.show_file() else "")

    def render_atoms_ro(self):
        self.atoms_ro_bridge.value = self.atoms_ro_html()

    # --- grid selection (per-section) -------------------------------------------------
    def selected_grid_id(self):
        """The selected grid id, validated against the live registry plus the Individual
        pseudo-grid (falls back to the first registered grid, else Individual)."""
        ids = [gid for gid, _ in self.grids.names] + [INDIVIDUAL_GRID]
        cur = store.get_section(self.grid_section).get("atmos_grid_id")
        if cur not in ids:
            cur = ids[0]
            store.update_pref(self.grid_section, "atmos_grid_id", cur)
        return cur

    @property
    def is_individual(self):
        return self.selected_grid_id() == INDIVIDUAL_GRID

    @property
    def is_nlte(self):
        return str(self.mode_form.value("mode")).upper() == "NLTE"

    def selected_grid(self):
        """The selected *registered* grid dict, or ``None`` (no grids, or Individual)."""
        return None if self.is_individual else self.grids._find(self.selected_grid_id())

    def active_grid(self):
        """The grid dict to synthesise with: a registered grid, or a synthetic dict for the
        Individual pseudo-grid (no folder/pattern; format from the Format selector)."""
        if self.is_individual:
            fmt = store.get_section(self.section).get("atmos_individual_format") or ATMOS_FORMATS[0]
            return {"name": "Individual", "folder": "", "pattern": "", "format": fmt}
        return self.grids._find(self.selected_grid_id())

    def atmos_file(self):
        """The concrete atmosphere file path when one is in play (Individual or a picked grid
        member), else ``None`` — the caller then interpolates the grid."""
        if not self.show_file():
            return None
        return (store.get_section(self.section).get("atmos_file") or "").strip() or None

    def _individual_dir(self):
        cwd = Path(str(load_prefs().get("dispatch_cwd") or REPO_DIR)).expanduser()
        d = cwd / "input" / "atmos"
        return d if d.is_dir() else cwd

    def _atmos_picker_root(self):
        if self.is_individual:
            return self._individual_dir()
        folder = ((self.selected_grid() or {}).get("folder") or "").strip()
        return Path(folder).expanduser() if folder else None

    def _reset_atmos_picker(self):
        self.atmos_picker.set_root(self._atmos_picker_root())

    # --- atmosphere UI rendering ------------------------------------------------------
    def grid_ui_html(self):
        """Grid selector (registered grids + 'Individual atmosphere'); then a Source toggle
        (Interpolate / Pick file) for a registered grid, or a Format selector for Individual.
        The chosen file itself lives in the always-on path field below this block (revealed by
        :meth:`show_file`), so no chip/browse button is rendered here. Namespaced by
        ``self.prefix``."""
        p = self.prefix
        sel = self.selected_grid_id()
        run = store.get_section(self.section)
        mode = run.get("atmos_mode", "interp")
        grid_opts = "".join(
            f'<option value="{html.escape(gid)}"{" selected" if gid == sel else ""}>'
            f'{html.escape(name)}</option>' for gid, name in self.grids.names)
        grid_opts += (f'<option value="{INDIVIDUAL_GRID}"'
                      f'{" selected" if sel == INDIVIDUAL_GRID else ""}>'
                      'Individual atmosphere</option>')
        out = ['<div class="field-group"><label>Grid</label>'
               f'<select class="field-input {p}-grid-select">{grid_opts}</select></div>']

        if sel == INDIVIDUAL_GRID:
            fmt = run.get("atmos_individual_format") or ATMOS_FORMATS[0]
            fmt_opts = "".join(
                f'<option value="{html.escape(f)}"{" selected" if f == fmt else ""}>'
                f'{html.escape(f)}</option>' for f in ATMOS_FORMATS)
            out.append('<div class="field-group field-inline"><label>Format</label>'
                       f'<select class="field-input {p}-format-select">{fmt_opts}</select></div>')
        else:
            mode_opts = "".join(
                f'<option value="{v}"{" selected" if v == mode else ""}>{lbl}</option>'
                for v, lbl in (("interp", "Interpolate"), ("pick", "Pick file")))
            out.append('<div class="field-group field-inline"><label>Source</label>'
                       f'<select class="field-input {p}-mode-select">{mode_opts}</select></div>')
        return "".join(out)

    def atoms_ro_html(self):
        """Read-only NLTE model-atom display (the selection is shared with the Synthesize
        tab, so it is managed there, not here)."""
        if not self.is_nlte:
            return ""
        atoms = list(self.atoms_provider())
        if not atoms:
            return ('<div class="field-static">NLTE selected — add model atoms in the '
                    'Synthesize tab; this run reuses them.</div>')
        chips = "".join(
            f'<span class="chip"><span class="chip-label">{html.escape(Path(a).name)}</span>'
            '</span>' for a in atoms)
        return ('<div class="hint">NLTE model atoms (shared with the Synthesize tab):</div>'
                f'<div class="chips">{chips}</div>')

    # --- action handling --------------------------------------------------------------
    def _on_grid_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "grid":
            store.update_pref(self.grid_section, "atmos_grid_id", arg)
            store.update_pref(self.section, "atmos_file", "")  # drop a stale cross-grid file
            self._reset_atmos_picker()
            self._seed_atmos_query()
            self.render_grid_ui()
        elif kind == "mode":
            store.update_pref(self.section, "atmos_mode", arg)
            self._seed_atmos_query()   # entering 'pick' shows the stored file (if any)
            self.render_grid_ui()
        elif kind == "format":
            store.update_pref(self.section, "atmos_individual_format", arg)
        elif kind == "clearfile":
            store.update_pref(self.section, "atmos_file", "")
            self._seed_atmos_query()
            self.render_grid_ui()

    def _on_pick_atmos(self, filepath):
        """A direct atmosphere file was picked: remember it and auto-fill the param fields
        from the filename via the grid's pattern (editable; non-matching files leave the
        fields as-is)."""
        path = str(filepath)
        store.update_pref(self.section, "atmos_file", path)
        self.atmos_picker.query.value = path   # the single field now shows the chosen file
        grid = self.selected_grid()
        if grid and grid.get("pattern"):
            rec = synth.recover_params(grid["pattern"], Path(path).name)
            if rec:
                teff, logg, feh, vmic = rec
                vals = {"teff": teff, "logg": logg, "feh": feh}
                if vmic is not None:
                    vals["vmic"] = vmic
                self.run_form.set_values(**vals)
                for k, v in vals.items():
                    store.update_pref(self.section, k, v)
        self.render_grid_ui()
        pulse(self.file_results_bridge, "off")   # collapse the dropdown, keep the field

    # --- persistence ------------------------------------------------------------------
    def _save_run(self, key, value):
        store.update_pref(self.section, key, value)

    def _save_base(self, key, value):
        store.update_pref(self.section, key, value)
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
                store.update_pref(self.section, "abundances", converted)
        if self.on_base_change:
            self.on_base_change(key, value)

    def _save_mode(self, key, value):
        store.update_pref(self.section, key, value)
        self.render_atoms_ro()
        if self.on_mode_change:
            self.on_mode_change(self.is_nlte)

    # --- inherit / reload -------------------------------------------------------------
    def inherit_from(self, run_section="synth_run", grid_section="synthesis",
                     grid_key="atmos_grid_id", *, force=False):
        """Seed this setup from another section's configuration. By default only fills keys
        this section does not already have (first-open inheritance); ``force=True`` overwrites
        (the explicit 'Copy from Synthesize tab' action)."""
        dest = dict(store.get_section(self.section))  # keeps the _seeded marker, if set
        src = store.get_section(run_section)
        for k in _INHERIT_KEYS:
            if k in src and (force or k not in dest):
                dest[k] = src[k]
        # always adopt the source's grid *selection* on inherit (the construction-time default
        # in this section must not suppress it); the Synthesize tab stores the id under
        # ``grid_section``/``grid_key`` ("synthesis"/"atmos_grid_id"), not its run section.
        gid = store.get_section(grid_section).get(grid_key)
        if gid:
            dest["atmos_grid_id"] = gid
        store.set_section(self.section, dest)
        self.reload()

    def reload(self):
        """Re-seed the forms + atmosphere UI from the (possibly just-changed) store section."""
        vals = {**RUN_DEFAULTS, **store.get_section(self.section)}
        self.run_form.set_values(**{k: vals[k] for k, _, _, _ in synth.STELLAR_PARAMS})
        self.base_form.set_values(abundances=vals["abundances"],
                                  relative_abundances=vals["relative_abundances"])
        self.mode_form.set_values(mode=vals["mode"])
        self._reset_atmos_picker()
        self._seed_atmos_query()
        self.render_grid_ui()
        self.render_atoms_ro()

    # --- the payload ------------------------------------------------------------------
    def context(self):
        """A :class:`~tumult.synthesis.SynthRequest` for the headless workflows, or ``None``
        (with :attr:`reason` set) when synthesis is not yet runnable."""
        prefs = load_prefs()
        if missing_prefs(prefs):
            self.reason = "Synthesis paths are missing — check Preferences → Synthesis."
            return None
        grid = self.active_grid()
        if grid is None:
            self.reason = ("No atmosphere grid — register one in Preferences → Synthesis, "
                           "or choose 'Individual atmosphere' in the Setup tab.")
            return None
        rv = self.run_form.values()
        rv.update(self.base_form.values())
        rv["mode"] = self.mode_form.value("mode")
        rv["atoms"] = list(self.atoms_provider())
        if self.is_nlte and not rv["atoms"]:
            self.reason = "NLTE needs model atoms — add them in the Synthesize tab."
            return None
        sect = store.get_section(self.section)
        atmos_file = None
        if self.is_individual or sect.get("atmos_mode", "interp") == "pick":
            atmos_file = (sect.get("atmos_file") or "").strip() or None
            if atmos_file is None:
                self.reason = ("Choose a grid to interpolate, or pick an atmosphere file, "
                               "in the Setup tab.")
                return None
        target = (float(rv["teff"]), float(rv["logg"]), float(rv["feh"]), float(rv["vmic"]))
        self.reason = None
        return SynthRequest(prefs=prefs, grid=grid, run_vals=rv,
                            target=target, atmos_file=atmos_file)
