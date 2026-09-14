"""The synthesis modal: one dialog that synthesises or fits a spectrum over a window.

It is opened on a *job* — a free wavelength range drawn with the Select range tool, the windows of
a line-mask set, or a blank window from the empty landing page — in one of two modes:

* **synth**: compute a spectrum at the configured (Teff, log g, [Fe/H], vmic), abundances
  and RT mode, and add it to the canvas. Every synthesis, whatever its atmosphere or
  window, is one durable workflow (``ctx.submit_workflow``): it is detached from the
  document, so closing the tab leaves it running and a later document can recover it.
* **fit**: fit a target spectrum by adaptive grid synthesis, starting from that same
  configuration, with the Fit section (:class:`fittab.FitTab`) saying what varies. The
  best-fit stellar parameters are written back into the setup, so the next synthesis starts
  from them. A fit is a durable ``fit`` workflow (``ctx.submit_task``) like a synthesis;
  ``_import_fit`` delivers a record another document submitted.

Both modes collect one owned :class:`~tumult.synthesis.NativeRequest` through
:meth:`synthsetup.SynthSetup.context` and add the job's windows to it: a free range is one
:class:`~tumult.synthesis.MaskWindow` whose fit and synth bounds coincide, a mask set
contributes its included masks with their distinct fit and synth windows. Unsupported
combinations (saved departures outside NLTE or with a grid, a fit with a snapshot
or saved departures, ...) are refused before anything is submitted, with the
reason in the footer and the dialog left open.

The atmosphere + stellar-parameter + RT-mode cluster is a composed
:class:`synthsetup.SynthSetup` (the same component the abundance / correction wizards own)
parameterised with this modal's historical store split (run state under ``"synth_run"``,
the grid selection under ``"synthesis"``) and form names (the page slots ``synth-fields`` /
``synthabund-fields`` / ``synthmode-fields``). The modal also owns the app-global pieces
that live in Preferences (the Atmosphere grids page, the NLTE model-atom selection,
the preference sheets and the run-cache listing). See ``docs/GUIWorkflow.md``.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

import store
from tumult.synthesis import FitRequest, GridTarget, MaskWindow, check_request
from tumult.synthesis import workflows as synth
from atompicker import AtomPicker
from controls import Field, Form
from fittab import FitTab
from gridmanager import GridManager
from filepicker import FilePicker
from htmlbridge import pulse, toggle
from synthsetup import SynthSetup
from synthprefs import (
    PATH_PICKERS,
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
    def __init__(self, ctx, add, on_warn=None):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self._add = add          # (spec, name, meta=None, conv=None) -> add to the canvas
        ctx.task_importers["fit"] = self._import_fit    # a fit record from another document
        self.on_inputs_change = lambda: None
        self.on_grids_change = lambda: None   # the grid registry changed (Preferences)
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

        self.pref_forms["synthengine"].groups = [
            ("Execution", ("synth_output_dir", "workflow_runner", "max_threads")),
            ("Spectral sampling", ("R", "linelist_crop_aa")),
            ("Atmosphere mesh (horizontal columns apply to 3D)", ("nx", "ny", "nz")),
            ("Fit settings", ("vmic_step", "abund_step", "fit_max_rounds"))]
        # the path preferences are picked, not typed: each has an inline file picker beside
        # its page's form (a Bokeh root, so it lives in the template rather than in the
        # re-rendered form) and its text field is hidden. Folders may be picked where the
        # preference accepts one.
        self.pref_pickers, self._pref_results, self._pref_warn = {}, {}, {}
        root = self._prefs_root(prefs)
        for key, dirs in PATH_PICKERS.items():
            form = next(f for f in self.pref_forms.values() if key in f.fields)
            form.fields[key].hidden = True
            picker = FilePicker(lambda path, key=key: self._on_pick_pref(key, path),
                                container_id=f"pref-{key}-results", pick_dirs=dirs,
                                root=None if key == "dispatch_cwd" else root, start_dir=root)
            picker.layout.name = f"pref_{key}_picker"
            picker.preset(str(prefs.get(key) or ""))
            self.pref_pickers[key] = picker
            self._pref_results[key] = ctx.bridges.cls(f"#pref-{key}-picker", "results-open")
            self._pref_warn[key] = ctx.bridges.cls(f"#pref-{key}-block", "field-warn")
        self._slurm_bridge = ctx.bridges.cls("#pref-slurm_script-block", "off")
        toggle(self._slurm_bridge, prefs.get("workflow_runner") != "slurm")

        # atmosphere grids (Preferences → Atmosphere grids); the setup's grid selector,
        # source-mode toggle and direct-file picker read the current selection from it.
        # Downloads go under <dispatch_cwd>/input/atmos unless a storage folder is set, so
        # the bundled MARCS grid there is registered from the start.
        cwd = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
        self.grids = GridManager(default_pattern=synth.DEFAULT_PATTERN,
                                 on_change=self._on_grids_changed,
                                 default_storage=cwd / "input" / "atmos", ctx=ctx)

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
    def _on_grids_changed(self):
        """GridManager edited the registry: refresh the selector, the model listing and
        the file warnings, and tell the app (the line-strength reference grid follows)."""
        self.setup.refresh_grids()
        self._refresh_warnings(push=True)
        self.on_grids_change()

    # --- opening ------------------------------------------------------------------------
    def open_range(self, xmin=None, xmax=None, mode="synth"):
        """Open on a free window (Å); ``None`` leaves the fields blank for the user."""
        self.job = ("range", None)
        self.setup.element = ""
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
        self.setup.element = mset.element
        self._open(mode)

    def _open(self, mode):
        self.draft.begin()
        self.mode = mode
        self.setup.allow_departures = mode != "fit"
        self.setup.allow_atom_only = mode != "fit"
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
            what = (f'<b>{len(masks)}</b> window{"" if len(masks) == 1 else "s"} of '
                    f'<b>{html.escape(mset.name)}</b>, {lo:.2f}–{hi:.2f} Å')
            if self.mode == "fit":
                return (f'<div class="hint">Fitting stellar parameters, abundances or '
                        f'broadening to a target spectrum, compared inside the fit windows '
                        f'of {what}.</div>')
            return (f'<div class="hint">Synthesising the synth windows of {what}, on the '
                    f'atmosphere chosen below.</div>')
        return ('<div class="section-header">Wavelength window</div>'
                '<div class="hint">Last-used atmosphere setup. Synthesis continues after closing the tab; spectrum fits stop.</div>'
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
                ["Microturbulence ξ [km/s]", f"{float(stellar['vmic']):.2f}"],
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
        """Fit the chosen target as a durable ``fit`` workflow (``ctx.submit_task``).

        The request is the setup's, with the job's windows; :class:`FitRequest` refuses
        what fitting does not support (a snapshot, saved departures, stellar axes
        with a fixed file) before the job is admitted."""
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason or "Cannot fit — check Preferences → Synthesis.")
            return
        windows = self._job_windows()
        if not windows:
            self.notify("Set the wavelength window (λ min/max) first.")
            return
        c = self.fit.collect()
        obj = self.canvas._find(c["target"]) if c["target"] else None
        if obj is None:
            self.notify("Select a target spectrum to fit.")
            return
        try:
            fit_request = FitRequest(request.with_masks(windows),
                                     fit_stellar=c["fit_stellar"],
                                     abund_fits=c["abund_fits"], fit_post=c["fit_post"])
        except ValueError as exc:
            self.notify(f"Cannot fit: {exc}")
            return
        synthesis = fit_request.synthesis
        rv, grid = synthesis.run_vals, self.setup.active_grid()
        guesses = {k: float(rv[k]) for k in ("teff", "logg", "feh", "vmic")}
        abundances, relative = rv["abundances"], bool(rv["relative_abundances"])
        observation = obj.spectrum
        xmin, xmax = fit_request.synth_span
        name = f"fit·{obj.name}"

        def on_done(result):
            best, raw, post = result
            self.fit.populate(best)
            # the best-fit stellar parameters become the setup, so the next synthesis (or
            # fit) starts from them; the broadening (post) is applied as a live convolution,
            # shown in the Edit tab.
            stellar = {k: float(best.get(k, guesses[k])) for k in guesses}
            self.adopt_stellar(stellar)
            terms = synth.fitted_abundance_terms(
                best, fit_request.abund_fits, base=abundances, relative=relative,
                feh=float(stellar["feh"]), solar=synth.abund_table(synthesis.prefs))
            meta = self._synth_meta(stellar, mode=rv["mode"], abundances=terms,
                                    relative=relative, atoms=rv["atoms"],
                                    xmin=xmin, xmax=xmax,
                                    atmosphere=self._selection_rows(synthesis.selection, grid))
            self._add(raw, name, meta=meta, conv=post)

        self.draft.commit()
        self.fit.commit()
        pulse(self.dialog, "close")
        # the fit only reads the observation inside the synth span (workflows.fit crops
        # it), so the record carries that part, with the guesses as the annotations a
        # document that imports the record later shows
        meta = self._synth_meta(guesses, mode=rv["mode"], abundances=abundances,
                                relative=relative, atoms=rv["atoms"], xmin=xmin, xmax=xmax,
                                atmosphere=self._selection_rows(synthesis.selection, grid))
        args = {"observation": observation.crop(xmin=xmin, xmax=xmax).to_dict(),
                "request": fit_request.to_dict()}
        self.ctx.submit_task(name, "fit", args, on_done, prefs=synthesis.prefs, meta=meta)

    def _import_fit(self, result, handle):
        """A fit record another document submitted: its best-fit spectrum is added under
        the recorded name and annotations, with the broadening as a live convolution."""
        best, raw, post = result
        self.fit.populate(best)
        self._add(raw, handle.label, meta=handle.meta, conv=post)

    def _job_windows(self):
        """The job's windows as :class:`MaskWindow` objects: a mask set's included masks
        with their own fit and synth windows (a mask without synth bounds spans its fit
        window), or the free range as one window with coinciding bounds. Empty when the
        range fields do not hold a valid interval."""
        kind, mset = self.job or ("range", None)
        if kind == "masks":
            windows = []
            for m in mset.masks:
                if not m.include:
                    continue
                synth_window = ((m.smin, m.smax) if m.smin == m.smin and m.smax == m.smax
                                else (m.lmin, m.lmax))
                windows.append(MaskWindow((m.lmin, m.lmax), synth_window))
            return tuple(windows)
        xmin, xmax = self._window_fields()
        return () if xmin is None else (MaskWindow((xmin, xmax), (xmin, xmax)),)

    def adopt_stellar(self, stellar):
        """Write ``{teff, logg, feh, vmic}`` into the setup (form + store)."""
        self.setup.set_stellar(stellar)
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
        prefs = load_prefs()
        if key == "dispatch_cwd":
            # the other path pickers browse relative to the working directory
            root = self._prefs_root(prefs)
            for other, picker in self.pref_pickers.items():
                if other != "dispatch_cwd":
                    picker.set_root(root)
                    picker.preset(str(prefs.get(other) or ""))
        if key in ("linelist", "dispatch_cwd"):
            self.on_inputs_change()
        if key == "max_threads":
            self._configure_daemon(prefs)
        toggle(self._slurm_bridge, prefs.get("workflow_runner") != "slurm")
        self._refresh_warnings(push=True)

    @staticmethod
    def _configure_daemon(prefs):
        """The Engine thread preference is the workflow daemon's budget; tell the daemon
        through this server's connection (a harness without one has no daemon to tell)."""
        from tumult.synthesis import daemon
        connection = daemon.current()
        if connection is None or connection.role != "server":
            return
        try:
            connection.request("configure", max_threads=int(prefs.get("max_threads") or 0))
        except (OSError, ValueError, daemon.DaemonError):
            logging.getLogger(__name__).warning("The daemon budget was not changed",
                                                exc_info=True)

    @staticmethod
    def _prefs_root(prefs):
        return Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()

    def _on_pick_pref(self, key, path):
        """A path preference was picked: show it in its field, persist it, fold the list."""
        value = str(path)
        self.pref_pickers[key].query.value = value
        self._save_pref(key, value)
        pulse(self._pref_results[key], "off")

    def pref_picker_results(self):
        """``{key: results-bridge name}`` for the template's inline pickers."""
        return {key: picker.results_name for key, picker in self.pref_pickers.items()}

    # --- file-existence warnings ------------------------------------------
    def _refresh_warnings(self, push):
        """Mark path fields whose file is missing; optionally re-render and notify."""
        self.missing = missing_prefs(load_prefs())
        for form in self.pref_forms.values():
            for key, field in form.fields.items():
                field.warn = key in self.missing
            if push:
                form.render()
        for key, bridge in self._pref_warn.items():
            toggle(bridge, key in self.missing)
        if push and self.on_warn:
            self.on_warn(self.has_warning)

    # --- launch -----------------------------------------------------------
    def synthesize(self):
        """Submit the job as one durable synthesis workflow, whatever its atmosphere."""
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason
                        or "Cannot synthesize — check Preferences → Synthesis.")
            return
        windows = self._job_windows()
        if not windows:
            self.notify("Set the synthesis wavelength window (λ min/max) first.")
            return
        try:
            request = check_request(request.with_masks(windows))
        except ValueError as exc:
            self.notify(f"Cannot synthesize: {exc}")
            return
        kind, mset = self.job or ("range", None)
        rv, target = request.run_vals, request.target
        name = (f"t{target[0]:.0f}g{target[1]:.1f}m{target[2]:+.1f}v{target[3]:.1f}"
                + ("·NLTE" if str(rv.get("mode")).upper() == "NLTE" else ""))
        shown_min = min(w.synth[0] for w in windows)
        shown_max = max(w.synth[1] for w in windows)
        atmosphere = self._selection_rows(request.selection, self.setup.active_grid())
        if request.departures is not None:
            atmosphere.append(["Saved departures", str(request.departures.folder)])
        meta = self._synth_meta(
            {"teff": target[0], "logg": target[1], "feh": target[2], "vmic": target[3]},
            mode=rv["mode"], abundances=rv["abundances"],
            relative=rv["relative_abundances"], atoms=rv["atoms"],
            xmin=shown_min, xmax=shown_max, atmosphere=atmosphere)
        linelist = ("Off (model atom only)" if rv.get("only_atom") else
                    str(request.prefs.get("linelist", ""))
                    + (f" ({rv['only_elements']} lines only)" if rv.get("only_elements") else ""))
        meta["params"].append(["Line list", linelist])
        if kind == "masks":
            meta["params"].append(["Line mask", mset.name])

        def add_result(result):
            for index, spec in enumerate(result.spectra, start=1):
                suffix = f"·{index}" if len(result.spectra) > 1 else ""
                self._add(spec, name + suffix, meta=meta)

        self.draft.commit()
        self.fit.commit()
        pulse(self.dialog, "close")
        self.ctx.submit_workflow(name, request, add_result, meta)
