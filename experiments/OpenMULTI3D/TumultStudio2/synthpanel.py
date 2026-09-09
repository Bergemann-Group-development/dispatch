"""The Synthesize tool: a toolbar button that swaps the right sidebar to a small
run-settings panel (Teff, logg, FeH, vmic, abundances, …) and computes a spectrum at
that point, plus the reusable synthesis preferences shown in the settings modal.

The atmosphere + stellar-parameter + RT-mode cluster is a composed
:class:`synthsetup.SynthSetup` (the same component the abundance / correction wizards
own) parameterised with this tab's historical store split (run state under
``"synth_run"``, the grid selection under ``"synthesis"``) and form names (the page
slots ``synth-fields`` / ``synthabund-fields`` / ``synthmode-fields``). This panel adds
what is unique to the Synthesize tab: the wavelength window, the Fit sub-tab, the run
cache and the preference sheets.
"""
from __future__ import annotations

import html
from pathlib import Path

import astropy.units as u

import store
import synth
from atompicker import AtomPicker
from controls import Field, Form
from fittab import FitTab
from gridmanager import GridManager
from rangeselect import RangeSelect
from synthsetup import SynthSetup
# shared prefs/constants live in synthprefs; re-exported here for back-compat with
# older imports (`from synthpanel import RUN_DEFAULTS, ...`).
from synthprefs import (  # noqa: F401
    INDIVIDUAL_GRID,
    INTERP_METHODS,
    PATH_PREFS,
    PREF_DEFAULTS,
    PREF_PAGES,
    REPO_DIR,
    RUN_DEFAULTS,
    SELECT_OPTIONS,
    load_prefs,
    missing_prefs,
)
from units import to_range

# the drawn synthesis band — green, to distinguish it from the EW tool's blue selection
SYNTH_BAND = "#50c878"


def _human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


class SynthPanel(RangeSelect):
    def __init__(self, canvas, submit, add, on_mode=None, on_warn=None, notify=None):
        super().__init__(canvas, color=SYNTH_BAND, description="Select synthesis x-range")
        self._submit = submit    # (name, loader, on_done) -> runs loader on the synth executor
        self._add = add          # (spec, name) -> add the result to the canvas
        self.on_mode = on_mode   # push "on"/"off" to show/hide the panel
        self.on_warn = on_warn   # push True/False when configured files go missing/found
        self.notify = notify or (lambda *a, **k: None)  # footer info message
        self.active = False
        self.missing = []        # path prefs currently pointing at non-existent files

        # run-cache listing (settings modal → Cache page)
        from htmlbridge import html_bridge
        self.cache_bridge = html_bridge("synth-cache", name="synth-cache-render")

        # synthesis wavelength window (Synthesize tab only): seeded from the drawn selection /
        # visible view when the panel opens, and required before Synthesize will run. Left blank
        # on the empty landing page (the empty canvas has no meaningful range to seed from).
        self.window_form = Form("synthwindow", [
            Field("xmin", "λ min [Å]", kind="number", value=None, min=0, inline=True,
                  required=True),
            Field("xmax", "λ max [Å]", kind="number", value=None, min=0, inline=True,
                  required=True),
        ])

        prefs = load_prefs()
        self.pref_forms = {
            name: Form(name, [
                Field(key, label, kind=kind, value=prefs.get(key, default),
                      options=SELECT_OPTIONS.get(key))
                for key, label, kind, default in fields
            ], on_change=self._save_pref)
            for name, fields in PREF_PAGES.items()
        }

        # registered atmosphere grids (Preferences → Synthesis); the setup's grid selector,
        # source-mode toggle and direct-file picker read the current selection from it. The
        # bundled MARCS grid under <dispatch_cwd>/input/atmos/MARCS is auto-registered on
        # first launch when present, so a fresh checkout works after fetching the input data.
        cwd = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
        self.grids = GridManager(default_pattern=synth.DEFAULT_PATTERN,
                                 on_change=self._on_grids_changed,
                                 start_dir=self._grid_start_dir(),
                                 builtin_dir=cwd / "input" / "atmos")

        # NLTE model-atom selector: opens in <dispatch_cwd>/input/atoms, shown only in NLTE.
        atoms_dir = cwd / "input" / "atoms"
        self.atoms = AtomPicker(start_dir=atoms_dir if atoms_dir.is_dir() else cwd)

        # the shared atmosphere + stellar-parameter + RT-mode cluster, with this tab's
        # historical store split and page-slot form names.
        self.setup = SynthSetup(
            "synth", "synth_run", self.grids,
            atoms_provider=lambda: list(self.atoms.selected),
            notify=self.notify,
            on_mode_change=self.atoms.show,          # reveal the atom selector only for NLTE
            on_base_change=self._on_base_change,     # keep the Fit tab's rows in sync
            on_atmos_change=self._sync_fit_stellar,  # hide stellar switches for a fixed atmosphere
            grid_section="synthesis",
            form_names={"run": "synth", "base": "synthabund", "mode": "synthmode"})

        self._refresh_warnings(push=False)  # seed field.warn so the initial html() is correct

        # Fit tab: shares RT-mode + atoms + abundances from the Synthesize tab; targets come
        # from the canvas. Seed its display flag from the shared relative-abundances toggle.
        self.fit = FitTab(on_run=self.run_fit)
        self.fit.set_catalog(prefs)
        self.fit.set_relative(bool(self.setup.base_form.value("relative_abundances")))
        self.fit.targets = [(o.uid, o.name) for o in self.canvas.objects]
        self._sync_fit_stellar()   # seed the stellar-fit switches for the current atmosphere

    # --- the composed setup, exposed under the panel's historical API -------------------
    @property
    def run_form(self):
        return self.setup.run_form

    @property
    def abund_form(self):
        return self.setup.base_form

    @property
    def mode_form(self):
        return self.setup.mode_form

    @property
    def atmos_picker(self):
        return self.setup.atmos_picker

    @property
    def is_nlte(self):
        return self.setup.is_nlte

    def selected_grid(self):
        return self.setup.selected_grid()

    def grid_ui_html(self):
        return self.setup.grid_ui_html()

    @property
    def models(self):
        return [*self.setup.models, *self.window_form.models,
                self.cache_bridge, *self.atoms.models, *self.fit.models,
                *self.grids.models,
                *(m for form in self.pref_forms.values() for m in form.models)]

    # --- registered grids ---------------------------------------------------------------
    def _grid_start_dir(self):
        """Where the grid-folder picker first opens: the selected grid's folder (falling back
        to the first registered grid, then home). Reads the store directly so it is safe to
        call before ``self.setup`` exists."""
        grids = store.get_section("synthesis").get("atmos_grids", [])
        if not grids:
            return Path("~")
        gid = store.get_section("synthesis").get("atmos_grid_id")
        g = next((x for x in grids if x.get("id") == gid), grids[0])
        folder = (g.get("folder") or "").strip()
        return Path(folder).expanduser() if folder else Path("~")

    def _on_grids_changed(self):
        """GridManager edited the registry: refresh the selector and file warnings."""
        self.setup.render_grid_ui()
        self._refresh_warnings(push=True)

    # --- Info-tab provenance ----------------------------------------------
    @staticmethod
    def _synth_meta(stellar, *, mode, abundances, relative, atoms, xmin, xmax):
        """A structured Info-tab payload describing the synthesis that produced a spectrum.

        ``stellar`` is ``{teff, logg, feh, vmic}``; the post-processing (vshift/vmac/vrot)
        is intentionally omitted here — it lives in the Edit tab as a live convolution.
        """
        rows = [["Teff [K]", f"{float(stellar['teff']):.0f}"],
                ["log g", f"{float(stellar['logg']):.2f}"],
                ["[Fe/H]", f"{float(stellar['feh']):+.2f}"],
                ["vmic [km/s]", f"{float(stellar['vmic']):.2f}"],
                ["RT mode", str(mode).upper()]]
        if str(abundances).strip():
            unit = "[X/Fe]" if relative else "A(X)"
            rows.append([f"Abundances ({unit})", str(abundances).strip()])
        if str(mode).upper() == "NLTE" and atoms:
            rows.append(["Model atoms", ", ".join(Path(a).name for a in atoms)])
        rows.append(["Window", f"{xmin:.2f}–{xmax:.2f} Å" if xmin and xmax else "full range"])
        # `stellar` is the same values as numbers rather than display strings: linestrength
        # asks an object which star it is, to pick the reference run its strengths come from
        return {"type": "synth", "params": rows,
                "stellar": {k: float(stellar[k]) for k in ("teff", "logg", "feh", "vmic")}}

    # --- fit tab ----------------------------------------------------------
    def refresh_targets(self):
        """Re-populate the Fit tab's target dropdown from the current canvas objects."""
        self.fit.set_targets([(o.uid, o.name) for o in self.canvas.objects])

    def run_fit(self, fittab):
        """Run an adaptive grid fit of the chosen target on the synth executor."""
        if self.missing:
            self.notify(f"Cannot fit — missing files: {', '.join(self.missing)}. Check Preferences.")
            return
        # A registered grid (interpolated) or a single fixed atmosphere file (Individual /
        # picked grid member). The fixed-file case fits abundances + broadening only — the Fit
        # tab already hides the stellar switches, so ``fit_stellar`` comes back empty.
        grid = self.setup.active_grid()
        if grid is None:
            self.notify("No atmosphere grid — register one in Preferences → Synthesis, or "
                        "choose 'Individual atmosphere' in the Setup tab.")
            return
        atmos_file = self.setup.atmos_file()
        if self.setup.show_file() and not atmos_file:
            self.notify("Pick an atmosphere file in the Setup tab to fit against it.")
            return
        c = fittab.collect()
        obj = self.canvas._find(c["target"]) if c["target"] else None
        if obj is None:
            self.notify("Select a target spectrum to fit.")
            return
        if not c["fit_stellar"] and not c["abund_fits"] and not c["fit_post"]:
            self.notify("Switch on at least one parameter (or add an abundance) to fit.")
            return
        mode = self.mode_form.value("mode")
        if str(mode).upper() == "NLTE" and not self.atoms.selected:
            self.notify("NLTE fit needs model atoms — add them via the RT mode selector.")
            return
        av = self.abund_form.values()  # shared Abundances + relative flag (one source of truth)
        abundances, relative = av["abundances"], bool(av["relative_abundances"])
        run_vals = {**c["guesses"], "abundances": abundances,
                    "relative_abundances": relative, "mode": mode,
                    "atoms": list(self.atoms.selected)}
        target = obj.spectrum
        prefs = load_prefs()
        xmin, xmax = self._window_aa()
        name = f"fit·{obj.name}"

        def on_done(result):
            best, raw, post = result
            fittab.populate(best)
            # Info-tab provenance: the best-fit stellar params + fitted abundances. The
            # broadening (post) is applied as a live convolution, shown in the Edit tab.
            stellar = {k: best.get(k, c["guesses"][k]) for k in ("teff", "logg", "feh", "vmic")}
            solar = synth.abund_table(prefs)
            feh = float(stellar["feh"])
            terms = [abundances.strip()] if str(abundances).strip() else []
            for el in c["abund_fits"]:
                if el in best:
                    xfe = float(best[el])
                    val = xfe if relative else synth.xfe_to_abund(xfe, el, feh, solar)
                    terms.append(f"{el}={val:+.2f}")
            meta = self._synth_meta(stellar, mode=mode, abundances=" ".join(terms),
                                    relative=relative, atoms=run_vals["atoms"],
                                    xmin=xmin, xmax=xmax)
            self._add(raw, name, meta=meta, conv=post)

        self._submit(
            name,
            lambda token=None: synth.fit(target, prefs, run_vals, c["fit_stellar"],
                                         c["abund_fits"], c["fit_post"], xmin, xmax,
                                         grid=grid, atmos_file=atmos_file, cancel_token=token),
            on_done)

    # --- run cache --------------------------------------------------------
    @staticmethod
    def _abund_text(abund, feh):
        """Per-run abundance overrides as '[El/Fe] +0.30' chips. Absolute A(X) values
        (relative_abundances off) are converted to [X/Fe] via the solar table + [Fe/H]."""
        items = abund.get("items") if abund else None
        if not items:
            return ""
        solar = {} if abund["relative"] else synth.abund_table(load_prefs())
        parts = []
        for el, val in items:
            if abund["relative"]:
                xfe = val
            else:
                a0 = solar.get(el)
                if a0 is None:
                    continue  # no solar reference -> can't express as [X/Fe]
                xfe = val - a0 - (feh or 0.0)
            parts.append(f"[{el}/Fe] {xfe:+.2f}")
        return " · ".join(parts)

    def cache_html(self):
        entries = synth.list_cache(load_prefs())
        total = sum(e["size"] for e in entries)
        head = (f'<div class="cache-summary">{len(entries)} cached run'
                f'{"" if len(entries) == 1 else "s"} · {_human_size(total)}</div>')
        if not entries:
            return head + '<div class="field-static">Cache is empty.</div>'
        rows = []
        for e in entries:
            star = e["star"]
            is_dep = e["kind"] == "departure"
            tag = ('<span class="cache-tag cache-tag-dep">DEP</span>' if is_dep
                   else '<span class="cache-tag">SPEC</span>')
            title = (f"{star['teff']:.0f} K · log g {star['logg']:.2f} · [Fe/H] {star['feh']:+.2f}"
                     + (f" · ξ {star['vmic']}" if star and star["vmic"] else "")) if star else e["run"]
            win = "departures" if is_dep else (
                f"{e['window'][0]}–{e['window'][1]} Å" if e["window"] else "full range")
            abund = self._abund_text(e["abund"], star["feh"] if star else 0.0)
            sub = " · ".join(filter(None, [win, abund, _human_size(e["size"])]))
            rows.append(
                f'<div class="cache-item">'
                f'<div class="cache-meta">'
                f'<span class="cache-star">{tag}{html.escape(title)}</span>'
                f'<span class="cache-sub">{html.escape(sub)}</span></div>'
                f'<button class="cache-del" data-path="{html.escape(e["path"])}" '
                f'title="Delete this cached run"><i class="ti ti-trash"></i></button>'
                f'</div>')
        return head + "".join(rows)

    def refresh_cache(self):
        self.cache_bridge.value = self.cache_html()

    def delete_cache(self, path):
        synth.delete_cache_entry(load_prefs(), path)
        self.refresh_cache()

    def clear_cache(self):
        synth.clear_cache(load_prefs())
        self.refresh_cache()

    @property
    def has_warning(self):
        return bool(self.missing) or bool(self.grids.missing_folders())

    # --- persistence ------------------------------------------------------
    def _on_base_change(self, key, value):
        """The shared Abundances / relative flag changed (persisted by the setup); keep the
        Fit tab's per-element row display in the same convention."""
        if key == "relative_abundances":
            self.fit.set_relative(bool(value))

    def _sync_fit_stellar(self):
        """Stellar params can only be fit by interpolating a grid; when the atmosphere is a
        concrete file (Individual, or a picked grid member) the Fit tab hides those switches."""
        if hasattr(self, "fit"):
            self.fit.set_can_fit_stellar(not self.setup.show_file())

    def _save_pref(self, key, value):
        store.update_pref("synthesis", key, value)
        self._refresh_warnings(push=True)

    # --- file-existence warnings ------------------------------------------
    def _refresh_warnings(self, push):
        """Mark path fields whose file is missing; optionally re-render and notify."""
        self.missing = missing_prefs(load_prefs())
        for form in self.pref_forms.values():
            for key, field in form.fields.items():
                field.warn = key in self.missing
            if push:
                form.render()
        if push and self.on_warn:
            self.on_warn(self.has_warning)

    # --- panel toggle + x-range selection ---------------------------------
    def toggle(self):
        """Synth button: open the panel and arm the x-range select tool, or close it."""
        self.cancel() if self.active else self.activate()

    def activate(self):
        self.active = True
        self.arm()                       # arm the box-select tool (RangeSelect)
        self._seed_window()              # fill λ min/max from the current view (blank if empty)
        if self.on_mode:
            self.on_mode("on")

    def _seed_window(self):
        """Fill the λ min/max fields from the drawn selection, else the visible view (both in
        Å). With nothing loaded the canvas only has a placeholder range, so leave them blank —
        the user must set a window before synthesising (see ``synthesize``)."""
        lo = hi = None
        if self.selection is not None:
            lo, hi = to_range(self.selection[0], self.selection[1], self.canvas.xunit, u.AA)
        elif self.canvas.objects:
            r = self.main.x_range
            if r.start is not None and r.end is not None:
                lo, hi = to_range(r.start, r.end, self.canvas.xunit, u.AA)
        if lo is not None and hi is not None:
            self.window_form.set_values(xmin=round(min(lo, hi), 3), xmax=round(max(lo, hi), 3))
        else:
            self.window_form.set_values(xmin=None, xmax=None)

    def on_select(self, x0, x1):
        """A range was drawn on the canvas: shade it and mirror it into the λ min/max fields."""
        self.show_band(x0, x1)
        lo, hi = to_range(self.selection[0], self.selection[1], self.canvas.xunit, u.AA)
        self.window_form.set_values(xmin=round(min(lo, hi), 3), xmax=round(max(lo, hi), 3))

    def cancel(self):
        """Leave synth mode: drop the drawn range, hide the band, release the tool."""
        was = self.active or self.selection is not None
        self.active = False
        self.clear_selection()
        self.disarm()                    # restore panning
        if was and self.on_mode:
            self.on_mode("off")

    # uniform ToolManager interface (the synth-cancel action / Esc path uses cancel())
    deactivate = cancel

    def _window_aa(self):
        """Synthesis window in Angstrom: the drawn range, else the current canvas view.
        Used by the Fit tab (the Synthesize tab reads its explicit fields instead)."""
        if self.selection is not None:
            x0, x1 = self.selection
        else:
            r = self.main.x_range
            if r.start is None or r.end is None:
                return None, None
            x0, x1 = r.start, r.end
        return to_range(x0, x1, self.canvas.xunit, u.AA)

    def _window_fields(self):
        """The Synthesize tab's λ min/max (Å), or (None, None) when either is blank or the
        pair is not a valid ascending, positive interval."""
        v = self.window_form.values()
        try:
            lo, hi = float(v["xmin"] or 0), float(v["xmax"] or 0)
        except (TypeError, ValueError):
            return None, None
        return (lo, hi) if hi > lo > 0 else (None, None)

    # --- launch -----------------------------------------------------------
    def synthesize(self):
        ctx = self.synth_context()
        if ctx is None:
            self.notify(self.setup.reason
                        or "Cannot synthesize — check Preferences → Synthesis.")
            return
        xmin, xmax = self._window_fields()
        if xmin is None:
            self.notify("Set the synthesis wavelength window (λ min/max) first.")
            return
        rv, grid, target = ctx.run_vals, ctx.grid, ctx.target
        atmos_file, prefs = ctx.atmos_file, ctx.prefs
        name = (f"t{target[0]:.0f}g{target[1]:.1f}m{target[2]:+.1f}v{target[3]:.1f}"
                + ("·NLTE" if str(rv.get("mode")).upper() == "NLTE" else ""))
        meta = self._synth_meta(
            {"teff": target[0], "logg": target[1], "feh": target[2], "vmic": target[3]},
            mode=rv["mode"], abundances=rv["abundances"],
            relative=rv["relative_abundances"], atoms=rv["atoms"], xmin=xmin, xmax=xmax)
        self._submit(name,
                     lambda token=None: synth.compute(target, prefs, rv, xmin, xmax,
                                                       grid=grid, atmos_file=atmos_file,
                                                       cancel_token=token),
                     lambda spec: self._add(spec, name, meta=meta))

    def synth_context(self):
        """Current synthesis configuration — delegated to the composed setup, which is the
        single source of truth (same grid, stellar params, RT mode, atoms and abundances).

        Returns ``None`` when synthesis is not runnable (``self.setup.reason`` says why);
        otherwise a dict ``{prefs, grid, run_vals, target, atmos_file}`` ready for
        ``synth.compute``.
        """
        return self.setup.context()
