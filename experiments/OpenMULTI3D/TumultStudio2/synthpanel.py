"""The Synthesize tool: a toolbar button that swaps the right sidebar to a small
run-settings panel (Teff, logg, FeH, vmic, abundances, …) and computes a spectrum at
that point, plus the reusable synthesis preferences shown in the settings modal.

Mirrors the EW tool's structure: a `controls.Form` for the panel, a class bridge to
show/hide it, and the heavy compute handed to a single-slot background executor.
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
from htmlbridge import html_bridge
from rangeselect import RangeSelect
from units import to_range

# the drawn synthesis band — green, to distinguish it from the EW tool's blue selection
SYNTH_BAND = "#50c878"

# the dispatch experiment directory (holds dispatch.x and input/); relative file
# prefs (the Fortran defaults below) are resolved against it for existence checks.
REPO_DIR = str(Path(synth.__file__).resolve().parent.parent)

# --- preferences: the reusable namelist setup, grouped into modal pages -----------------
# each page is one controls.Form; (key, label, kind, default) per field.
ATMOS_FORMATS = ["Marcs", "Stagger", "Stagger2", "Co5bold", "Text"]
INTERP_METHODS = ["NaturalNeighbor", "Delaunay"]   # NN: smoother (more nodes); Delaunay: fast
SELECT_OPTIONS = {"atmos_format": ATMOS_FORMATS, "interp_method": INTERP_METHODS}
# defaults follow the Fortran namelist defaults (composition/atom/atmos/spectrum _params);
# file paths may be relative to the Dispatch run directory (resolved for existence checks).
PREF_PAGES = {
    "synthprefs": [   # "Synthesis" page — the science inputs
        ("atmos_grid_folder", "Atmosphere grid folder", "text", "/Users/hoppe/Documents/MARCS"),
        ("atmos_grid_pattern", "Filename pattern", "text", synth.DEFAULT_PATTERN),
        ("atmos_format", "Atmosphere format", "select", "Marcs"),
        ("interp_method", "Interpolation", "select", "NaturalNeighbor"),
        ("linelist", "Synthesis line list (file or folder)", "text", ""),
        ("abund_file", "Abundance file", "text", "./input/abund/abund_magg"),
        ("absmet", "Absmet file or folder", "text", ""),
    ],
    "synthengine": [  # "Engine" page — paths + runtime knobs
        ("default_nml", "Default namelist", "text", ""),
        ("dispatch_cwd", "Dispatch run directory", "text", REPO_DIR),
        ("synth_output_dir", "Output directory", "text", "synth"),
        ("R", "Resolving power R", "number", 300000),
        ("nz", "Atmosphere depth points nz", "number", 128),
        ("vmic_step", "vmic lattice step [km/s]", "number", 0.5),
        # fit controls: Teff/logg/[Fe/H] are sampled at the atmosphere-grid nodes (no step);
        # only the continuous overrides (vmic via vmic_step, abundances) have a lattice step.
        ("abund_step", "Fit: abundance grid step [dex]", "number", 0.25),
        ("fit_max_rounds", "Fit: max refinement rounds", "number", 6),
        ("max_threads", "Max threads (0 = all cores)", "number", 0),
    ],
}
PREF_DEFAULTS = {key: default
                 for fields in PREF_PAGES.values()
                 for key, _, _, default in fields}

# path prefs that must exist on disk if non-empty (output dir is created on demand);
# relative ones resolve against dispatch_cwd.
PATH_PREFS = ["atmos_grid_folder", "linelist", "abund_file",
              "absmet", "default_nml", "dispatch_cwd"]


def _resolve_path(prefs, value):
    p = Path(str(value)).expanduser()
    if not p.is_absolute() and prefs.get("dispatch_cwd"):
        p = Path(str(prefs["dispatch_cwd"])).expanduser() / p
    return p


def missing_prefs(prefs):
    """Configured path prefs whose target does not exist (empty fields are ignored)."""
    missing = []
    for key in PATH_PREFS:
        val = prefs.get(key)
        if isinstance(val, str):
            val = val.strip()
        if not val:
            continue
        if not _resolve_path(prefs, val).exists():
            missing.append(key)
    return missing

# --- run settings: only the per-target knobs ----------------------------------------
RUN_DEFAULTS = {"teff": 5777, "logg": 4.4, "feh": 0.0, "vmic": 1.0,
                "abundances": "", "relative_abundances": True, "mode": "LTE",
                "atoms": []}


def load_prefs():
    """The persisted synthesis preferences merged over the defaults."""
    prefs = dict(PREF_DEFAULTS)
    prefs.update(store.get_section("synthesis"))
    return prefs


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
        self.cache_bridge = html_bridge("synth-cache", name="synth-cache-render")

        run_vals = dict(RUN_DEFAULTS)
        run_vals.update(store.get_section("synth_run"))
        self.run_form = Form("synth", [
            Field("teff", "Teff [K]", kind="number", value=run_vals["teff"], step=1, inline=True),
            Field("logg", "log g", kind="number", value=run_vals["logg"], step=0.01, inline=True),
            Field("feh", "[Fe/H]", kind="number", value=run_vals["feh"], step=0.01, inline=True),
            Field("vmic", "vmic [km/s]", kind="number", value=run_vals["vmic"], min=0, step=0.1,
                  inline=True),
            Field("abundances", "Abundances", kind="text", value=run_vals["abundances"],
                  placeholder="e.g. Ba=2.3, Eu=0.5"),
            Field("relative_abundances", "Relative abundances", kind="checkbox",
                  value=run_vals["relative_abundances"]),
        ], on_change=self._save_run)

        # RT mode lives in its own form so it can be rendered in a shared area used by both
        # the Synthesize and Fit tabs (one source of truth for LTE/NLTE + the atom selector).
        self.mode_form = Form("synthmode", [
            Field("mode", "RT mode", kind="select", value=run_vals["mode"],
                  options=["LTE", "NLTE"], inline=True),
        ], on_change=self._save_mode)

        prefs = load_prefs()
        self.pref_forms = {
            name: Form(name, [
                Field(key, label, kind=kind, value=prefs.get(key, default),
                      options=SELECT_OPTIONS.get(key))
                for key, label, kind, default in fields
            ], on_change=self._save_pref)
            for name, fields in PREF_PAGES.items()
        }
        self._refresh_warnings(push=False)  # seed field.warn so the initial html() is correct

        # NLTE model-atom selector: opens in <dispatch_cwd>/input/atoms, shown only in NLTE.
        cwd = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
        atoms_dir = cwd / "input" / "atoms"
        self.atoms = AtomPicker(start_dir=atoms_dir if atoms_dir.is_dir() else cwd)

        # Fit tab: shares RT-mode + atoms from the Synthesize tab; targets come from the canvas.
        self.fit = FitTab(on_run=self.run_fit)
        self.fit.set_catalog(prefs)
        self.fit.targets = [(o.uid, o.name) for o in self.canvas.objects]

    @property
    def is_nlte(self):
        return str(self.mode_form.value("mode")).upper() == "NLTE"

    @property
    def models(self):
        return [*self.run_form.models, *self.mode_form.models, self.cache_bridge,
                *self.atoms.models, *self.fit.models,
                *(m for form in self.pref_forms.values() for m in form.models)]

    # --- fit tab ----------------------------------------------------------
    def refresh_targets(self):
        """Re-populate the Fit tab's target dropdown from the current canvas objects."""
        self.fit.set_targets([(o.uid, o.name) for o in self.canvas.objects])

    def run_fit(self, fittab):
        """Run an adaptive grid fit of the chosen target on the synth executor."""
        if self.missing:
            self.notify(f"Cannot fit — missing files: {', '.join(self.missing)}. Check Preferences.")
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
        run_vals = {**c["guesses"], "abundances": c["abundances"],
                    "relative_abundances": c["relative"], "mode": mode,
                    "atoms": list(self.atoms.selected)}
        target = obj.spectrum
        prefs = load_prefs()
        xmin, xmax = self._window_aa()
        name = f"fit·{obj.name}"
        self._submit(
            name,
            lambda: synth.fit(target, prefs, run_vals, c["fit_stellar"],
                              c["abund_fits"], c["fit_post"], xmin, xmax),
            lambda result: (fittab.populate(result[0]), self._add(result[1], name)))

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
        return bool(self.missing)

    # --- persistence ------------------------------------------------------
    def _save_run(self, key, value):
        store.update_pref("synth_run", key, value)
        if key == "relative_abundances":
            # convert the fixed abundance overrides so the physical abundance is preserved
            text = str(self.run_form.value("abundances") or "").strip()
            if text:
                feh = float(self.run_form.value("feh") or 0.0)
                solar = synth.abund_table(load_prefs())
                converted = synth.convert_abundances(text, to_relative=bool(value),
                                                     feh=feh, solar=solar)
                self.run_form.set_values(abundances=converted)
                store.update_pref("synth_run", "abundances", converted)

    def _save_mode(self, key, value):
        store.update_pref("synth_run", key, value)
        self.atoms.show(self.is_nlte)  # reveal the model-atom selector only for NLTE

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
        if self.on_mode:
            self.on_mode("on")

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
        """Synthesis window in Angstrom: the drawn range, else the current canvas view."""
        if self.selection is not None:
            x0, x1 = self.selection
        else:
            r = self.main.x_range
            if r.start is None or r.end is None:
                return None, None
            x0, x1 = r.start, r.end
        return to_range(x0, x1, self.canvas.xunit, u.AA)

    # --- launch -----------------------------------------------------------
    def synthesize(self):
        if self.missing:
            self.notify("Cannot synthesize — missing files: "
                        f"{', '.join(self.missing)}. Check Preferences.")
            return
        rv = self.run_form.values()
        rv["mode"] = self.mode_form.value("mode")  # RT mode lives in its own shared form
        rv["atoms"] = list(self.atoms.selected)  # the NLTE model-atom selection (the form has no field for it)
        if self.is_nlte and not rv["atoms"]:
            self.notify("Cannot synthesize NLTE — no model atoms selected. Add one in the Synth panel.")
            return
        target = (float(rv["teff"]), float(rv["logg"]), float(rv["feh"]), float(rv["vmic"]))
        xmin, xmax = self._window_aa()
        prefs = load_prefs()
        name = (f"t{target[0]:.0f}g{target[1]:.1f}m{target[2]:+.1f}v{target[3]:.1f}"
                + ("·NLTE" if str(rv.get("mode")).upper() == "NLTE" else ""))
        self._submit(name,
                     lambda: synth.compute(target, prefs, rv, xmin, xmax),
                     lambda spec: self._add(spec, name))
