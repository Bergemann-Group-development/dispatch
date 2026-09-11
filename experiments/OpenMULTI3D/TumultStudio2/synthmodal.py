"""The synthesis modal: one dialog that synthesises or fits a spectrum over a window.

It is opened on a *job* — a free wavelength range drawn with the Mask tool, the windows of
a line-mask set, or a blank window from the empty landing page — in one of two modes:

* **synth**: compute a spectrum at the configured (Teff, log g, [Fe/H], vmic), abundances
  and RT mode, and add it to the canvas.
* **fit**: fit a target spectrum by adaptive grid synthesis, starting from that same
  configuration, with the Fit section (:class:`fittab.FitTab`) saying what varies. The
  best-fit stellar parameters are written back into the setup, so the next synthesis starts
  from them.

The atmosphere + stellar-parameter + RT-mode cluster is a composed
:class:`synthsetup.SynthSetup` (the same component the abundance / correction wizards own)
parameterised with this modal's historical store split (run state under ``"synth_run"``,
the grid selection under ``"synthesis"``) and form names (the page slots ``synth-fields`` /
``synthabund-fields`` / ``synthmode-fields``). The modal also owns the app-global pieces
that live in Preferences → Synthesis: the grid registry, the NLTE model-atom selection,
the preference sheets and the run-cache listing. See ``docs/GUIWorkflow.md``.
"""
from __future__ import annotations

import html
from pathlib import Path

import store
from tumult.synthesis import GridTarget, Snapshot3D, submit as submit_synthesis
from tumult.synthesis import workflows as synth
from atompicker import AtomPicker
from controls import Field, Form
from fittab import FitTab
from gridmanager import GridManager
from htmlbridge import pulse
from synthsetup import SynthSetup
from synthprefs import (
    PREF_PAGES,
    REPO_DIR,
    SELECT_OPTIONS,
    load_prefs,
    missing_prefs,
)


def _human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


class SynthModal:
    def __init__(self, ctx, submit, add, on_warn=None):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self._submit = submit    # (name, loader, on_done, detached=False) -> background run
        self._add = add          # (spec, name, meta=None, conv=None) -> add to the canvas
        self.on_warn = on_warn   # push True/False when configured files go missing/found
        self.missing = []        # path prefs currently pointing at non-existent files
        self.job = None          # ("range", None) | ("masks", LineMaskSet) once opened
        self.mode = "synth"      # "synth" | "fit"

        self.dialog = ctx.bridges.dialog("synth-modal", name="synth-modal-bridge")
        self.title_bridge = ctx.bridges.html("synth-modal-title")
        self.context_bridge = ctx.bridges.html("synth-context")
        self.actions_bridge = ctx.bridges.html("synth-modal-actions")
        self.fit_mode_bridge = ctx.bridges.cls("#synth-modal", "fit-mode")
        # run-cache listing (settings modal → Cache page)
        self.cache_bridge = ctx.bridges.html("synth-cache", name="synth-cache-render")

        # the wavelength window of a range job (Å); seeded from the range that opened the
        # modal, blank from the landing page, and required before a run starts
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

        # the shared atmosphere + stellar-parameter + RT-mode cluster, with this modal's
        # historical store split and page-slot form names.
        self.setup = SynthSetup(
            "synth", "synth_run", self.grids,
            atoms_provider=lambda: list(self.atoms.selected),
            notify=self.notify,
            on_mode_change=self.atoms.show,          # reveal the atom selector only for NLTE
            on_base_change=self._on_base_change,     # keep the Fit section's rows in sync
            on_run_change=self._on_run_change,       # ... and their absolute display
            on_atmos_change=self._sync_fit_stellar,  # hide stellar switches for a fixed atmosphere
            grid_section="synthesis",
            form_names={"run": "synth", "base": "synthabund", "mode": "synthmode"})

        self._refresh_warnings(push=False)  # seed field.warn so the initial html() is correct

        # the Fit section: target + what varies; the guess is the setup itself
        self.fit = FitTab()
        self.fit.set_catalog(prefs)
        self.fit.set_relative(bool(self.setup.base_form.value("relative_abundances")))
        self.fit.set_feh(self.setup.run_form.value("feh"))
        self.fit.targets = [(o.uid, o.name) for o in self.canvas.objects]
        self._sync_fit_stellar()   # seed the stellar-fit switches for the current atmosphere

    @property
    def notify(self):
        return self.ctx.notify

    # --- the composed setup, exposed under the modal's API ------------------------------
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

    def grid_ui_html(self):
        return self.setup.grid_ui_html()

    @property
    def models(self):
        return [*self.setup.models, *self.window_form.models,
                *self.atoms.models, *self.fit.models, *self.grids.models,
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

    # --- opening ------------------------------------------------------------------------
    def open_range(self, xmin=None, xmax=None, mode="synth"):
        """Open on a free window (Å); ``None`` leaves the fields blank for the user."""
        self.job = ("range", None)
        if xmin is not None and xmax is not None:
            self.window_form.set_values(xmin=round(float(xmin), 3), xmax=round(float(xmax), 3))
        else:
            self.window_form.set_values(xmin=None, xmax=None)
        self._open(mode)

    def open_masks(self, mset, mode="synth"):
        """Open on a line-mask set's included windows."""
        if mset is None or not any(m.include for m in mset.masks):
            self.notify("This line mask has no switched-on windows.")
            return
        self.job = ("masks", mset)
        self._open(mode)

    def _open(self, mode):
        self.mode = mode
        self.setup.render()
        self.atoms.show(self.is_nlte)
        self.refresh_targets()
        self.fit.render()
        pulse(self.fit_mode_bridge, "on" if mode == "fit" else "off")
        self.title_bridge.value = "Fit spectrum" if mode == "fit" else "Synthesise spectrum"
        self.context_bridge.value = self.context_html()
        self.actions_bridge.value = self.actions_html()
        pulse(self.dialog, "open")

    def context_html(self):
        """The modal's header block: the window the run covers."""
        kind, mset = self.job or ("range", None)
        if kind == "masks":
            masks = [m for m in mset.masks if m.include]
            lo, hi = self._mask_span(masks)
            verb = "Fitting inside" if self.mode == "fit" else "Synthesising"
            return (f'<div class="hint">{verb} <b>{len(masks)}</b> window'
                    f'{"" if len(masks) == 1 else "s"} of <b>{html.escape(mset.name)}</b>, '
                    f'{lo:.2f}–{hi:.2f} Å.</div>')
        return ('<div class="section-header">Wavelength window</div>'
                f'<div id="{self.window_form.slot_id}" class="synth-window">'
                f'{self.window_form.html()}</div>')

    def actions_html(self):
        if self.mode == "fit":
            return '<button class="primary-btn synth-submit">Fit</button>'
        return '<button class="primary-btn synth-submit">Synthesise</button>'

    def submit(self):
        """The footer button: what the modal was opened for."""
        if self.mode == "fit":
            self.run_fit()
        else:
            self.synthesize()

    # --- Info-tab provenance ----------------------------------------------
    @staticmethod
    def _atmosphere_rows(grid, atmos_file=None, resolution=None):
        """Info-tab rows naming the model atmosphere a synthesis read: the format the
        Fortran parsed it with, and either the interpolated grid or the one file."""
        grid = grid or {}
        rows = []
        if atmos_file:
            rows.append(["Atmosphere", Path(str(atmos_file)).name])
            if grid.get("folder"):
                rows.append(["Atmosphere grid", f"{grid.get('name') or grid.get('id')} (picked file)"])
        elif grid:
            rows.append(["Atmosphere", f"{grid.get('name') or grid.get('id')} grid, interpolated"])
        rows.append(["Atmosphere format", str(grid.get("format") or "unknown")])
        if resolution:
            rows.append(["Mesh", "×".join(str(int(v)) for v in resolution)])
        return rows

    @classmethod
    def _selection_rows(cls, selection, grid):
        """:meth:`_atmosphere_rows` for a request's atmosphere selection; ``grid`` is the
        setup's active grid record, which names the grid a picked file came from."""
        if isinstance(selection, GridTarget):
            return cls._atmosphere_rows(selection.grid)
        return cls._atmosphere_rows({**(grid or {}), "format": selection.format},
                                    selection.path, getattr(selection, "resolution", None))

    @staticmethod
    def _synth_meta(stellar, *, mode, abundances, relative, atoms, xmin, xmax,
                    atmosphere=()):
        """A structured Info-tab payload describing the synthesis that produced a spectrum.

        ``stellar`` is ``{teff, logg, feh, vmic}``; ``atmosphere`` the rows from
        :meth:`_atmosphere_rows`. The post-processing (vshift/vmac/vrot) is intentionally
        omitted here — it lives in the Edit tab as a live convolution.
        """
        rows = [["Teff [K]", f"{float(stellar['teff']):.0f}"],
                ["log g", f"{float(stellar['logg']):.2f}"],
                ["[Fe/H]", f"{float(stellar['feh']):+.2f}"],
                ["vmic [km/s]", f"{float(stellar['vmic']):.2f}"],
                *atmosphere,
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

    # --- the window --------------------------------------------------------
    @staticmethod
    def _mask_span(masks):
        """The synthesis span covering ``masks``: their synth windows, falling back to the
        fit window where a mask carries no synth bounds."""
        los = [m.smin if m.smin == m.smin else m.lmin for m in masks]
        his = [m.smax if m.smax == m.smax else m.lmax for m in masks]
        return min(los), max(his)

    def _window_fields(self):
        """The λ min/max fields (Å), or (None, None) when either is blank or the pair is not
        a valid ascending, positive interval."""
        v = self.window_form.values()
        try:
            lo, hi = float(v["xmin"] or 0), float(v["xmax"] or 0)
        except (TypeError, ValueError):
            return None, None
        return (lo, hi) if hi > lo > 0 else (None, None)

    # --- fit -------------------------------------------------------------
    def refresh_targets(self):
        """Re-populate the Fit section's target dropdown from the current canvas objects."""
        self.fit.set_targets([(o.uid, o.name) for o in self.canvas.objects])

    def run_fit(self):
        """Run an adaptive grid fit of the chosen target on the synth executor."""
        if self.missing:
            self.notify(f"Cannot fit — missing files: {', '.join(self.missing)}. Check Preferences.")
            return
        # A registered grid (interpolated) or a single fixed atmosphere file (Individual /
        # picked grid member). The fixed-file case fits abundances + broadening only — the Fit
        # section already hides the stellar switches, so ``fit_stellar`` comes back empty.
        grid = self.setup.active_grid()
        if grid is None:
            self.notify("No atmosphere grid — register one in Preferences → Synthesis, or "
                        "choose 'Individual atmosphere'.")
            return
        if self.setup.is_snapshot:
            self.notify("Fitting is available for 1D atmospheres and grids only.")
            return
        atmos_file = self.setup.atmos_file()
        if self.setup.show_file() and not atmos_file:
            self.notify("Pick an atmosphere file to fit against it.")
            return
        c = self.fit.collect()
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
        xmin, xmax, windows = self._fit_windows()
        if xmin is None:
            self.notify("Set the wavelength window (λ min/max) first.")
            return
        guesses = {k: float(v) for k, v in self.run_form.values().items()}
        av = self.abund_form.values()  # shared Abundances + relative flag (one source of truth)
        abundances, relative = av["abundances"], bool(av["relative_abundances"])
        run_vals = {**guesses, "abundances": abundances,
                    "relative_abundances": relative, "mode": mode,
                    "atoms": list(self.atoms.selected)}
        target = obj.spectrum
        prefs = load_prefs()
        name = f"fit·{obj.name}"

        def on_done(result):
            best, raw, post = result
            self.fit.populate(best)
            # the best-fit stellar parameters become the setup, so the next synthesis (or
            # fit) starts from them; the broadening (post) is applied as a live convolution,
            # shown in the Edit tab.
            stellar = {k: float(best.get(k, guesses[k])) for k in ("teff", "logg", "feh", "vmic")}
            self.adopt_stellar(stellar)
            terms = synth.fitted_abundance_terms(
                best, c["abund_fits"], base=abundances, relative=relative,
                feh=float(stellar["feh"]), solar=synth.abund_table(prefs))
            meta = self._synth_meta(stellar, mode=mode, abundances=terms,
                                    relative=relative, atoms=run_vals["atoms"],
                                    xmin=xmin, xmax=xmax,
                                    atmosphere=self._atmosphere_rows(grid, atmos_file))
            self._add(raw, name, meta=meta, conv=post)

        pulse(self.dialog, "close")
        self._submit(
            name,
            lambda token=None: synth.fit(target, prefs, run_vals, c["fit_stellar"],
                                         c["abund_fits"], c["fit_post"], xmin, xmax,
                                         windows=windows, grid=grid, atmos_file=atmos_file,
                                         cancel_token=token),
            on_done)

    def _fit_windows(self):
        """``(xmin, xmax, windows)`` for the fit: the synthesis span, and the fit windows the
        target is compared inside (``None`` for a free range, which is compared throughout)."""
        kind, mset = self.job or ("range", None)
        if kind == "masks":
            masks = [m for m in mset.masks if m.include]
            xmin, xmax = self._mask_span(masks)
            return xmin, xmax, [(m.lmin, m.lmax) for m in masks]
        xmin, xmax = self._window_fields()
        return xmin, xmax, None

    def adopt_stellar(self, stellar):
        """Write ``{teff, logg, feh, vmic}`` into the setup (form + store)."""
        self.run_form.set_values(**stellar)
        for k, v in stellar.items():
            store.update_pref("synth_run", k, v)
        self.fit.set_feh(stellar["feh"])

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
            elif el not in solar:
                continue  # no solar reference -> can't express as [X/Fe]
            else:
                xfe = synth.abund_to_xfe(val, el, feh or 0.0, solar)
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
        Fit section's per-element row display in the same convention."""
        if key == "relative_abundances":
            self.fit.set_relative(bool(value))

    def _on_run_change(self, key, value):
        if key == "feh":
            self.fit.set_feh(value)

    def _sync_fit_stellar(self):
        """Stellar params can only be fit by interpolating a grid; when the atmosphere is a
        concrete file (Individual, or a picked grid member) the Fit section hides those
        switches."""
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

    # --- launch -----------------------------------------------------------
    def synthesize(self):
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason
                        or "Cannot synthesize — check Preferences → Synthesis.")
            return
        kind, mset = self.job or ("range", None)
        selected_masks = tuple(m for m in mset.masks if m.include) if kind == "masks" else ()
        xmin, xmax = (None, None) if selected_masks else self._window_fields()
        if not selected_masks and xmin is None:
            self.notify("Set the synthesis wavelength window (λ min/max) first.")
            return
        rv, grid, target = request.run_vals, request.grid, request.target
        atmos_file, prefs = request.atmos_file, request.prefs
        name = (f"t{target[0]:.0f}g{target[1]:.1f}m{target[2]:+.1f}v{target[3]:.1f}"
                + ("·NLTE" if str(rv.get("mode")).upper() == "NLTE" else ""))
        shown_min, shown_max = (self._mask_span(selected_masks) if selected_masks
                                else (xmin, xmax))
        meta = self._synth_meta(
            {"teff": target[0], "logg": target[1], "feh": target[2], "vmic": target[3]},
            mode=rv["mode"], abundances=rv["abundances"],
            relative=rv["relative_abundances"], atoms=rv["atoms"],
            xmin=shown_min, xmax=shown_max,
            atmosphere=self._selection_rows(request.selection, self.setup.active_grid()))
        if selected_masks:
            meta["params"].append(["Line mask", mset.name])
        pulse(self.dialog, "close")
        if selected_masks or isinstance(request.selection, Snapshot3D):
            request_masks = selected_masks or ({"fit": (xmin, xmax), "synth": (xmin, xmax)},)
            request = request.with_masks(request_masks)

            def add_result(result):
                for index, spec in enumerate(result.spectra, start=1):
                    suffix = f"·{index}" if len(result.spectra) > 1 else ""
                    self._add(spec, name + suffix, meta=meta)

            def run_request(token=None):
                handle = submit_synthesis(request)
                detach = token.add_callback(handle.cancel) if token is not None else lambda: None
                try:
                    return handle.result()
                finally:
                    detach()

            self._submit(name, run_request, add_result, detached=True)
        else:
            self._submit(name,
                         lambda token=None: synth.compute(target, prefs, rv, xmin, xmax,
                                                           grid=grid, atmos_file=atmos_file,
                                                           cancel_token=token),
                         lambda spec: self._add(spec, name, meta=meta))
