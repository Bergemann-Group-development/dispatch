"""The synthesis modal: one dialog that synthesises or fits a spectrum over a window.

It is opened on a *job* — a free wavelength range drawn with the Select range tool, the windows of
a line-mask set, or a blank window from the empty landing page — in one of two modes:

* **synth**: compute a spectrum at the configured (Teff, log g, [Fe/H], vmic), abundances
  and RT mode, and add it to the canvas. Every synthesis, whatever its atmosphere or
  window, is one durable ``synthesis.synthesize`` workflow (``ctx.submit_operation``): it
  is detached from the document, so closing the tab leaves it running and a later document
  can recover it.
* **fit**: fit a target spectrum by adaptive grid synthesis, starting from that same
  configuration, with the Fit section (:class:`fittab.FitTab`) saying what varies. The
  best-fit stellar parameters are written back into the setup, so the next synthesis starts
  from them. A fit is a durable ``synthesis.fit`` workflow like a synthesis;
  ``_import_fit`` delivers a record another document submitted.

Both modes collect one :class:`~tumult.synthesis.SynthRequest` through
:meth:`synthsetup.SynthSetup.context` and give it the job's windows: a free range is one
window, a mask set contributes its included masks, whose synth windows are synthesised and
whose fit windows a fit compares inside. Unsupported combinations (saved departures with a
grid or for another star, a fit with a snapshot or saved departures, ...) are refused
before anything is submitted, with the reason in the footer and the dialog left open.

The atmosphere + stellar-parameter + RT-mode cluster is a composed
:class:`synthsetup.SynthSetup` (the same component the abundance / correction wizards own)
parameterised with this modal's historical store split (run state under ``"synth_run"``,
the grid selection under ``"synthesis"``) and form names (the page slots ``synth-fields`` /
``synthabund-fields`` / ``synthmode-fields``). The app-global pieces it shows but does not
own — the grid registry (Preferences → Atmosphere grids) and the NLTE model-atom selection
every calculation shares — are constructed before this dialog and passed in; the
preference sheets and the run-cache listing belong to :class:`preferences.page.PreferencesPage`.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

from ..workspace import store
from tumult.departures import Departures
from tumult.linemask import LineMask, LineMaskSet
from tumult.synthesis import AtmosphereGrid, Snapshot, SynthRequest, check_fit
from tumult.synthesis import workflows as synth
from ..ui.controls import Field, Form
from .fittab import FitTab
from ..ui.htmlbridge import pulse
from ..workspace import naming
from .synthsetup import SynthSetup
from ..preferences.synthprefs import load_prefs, synth_config


def rerun_record(request, windows, **extra):
    """What a spectrum's Info meta keeps to be rerun (``meta["rerun"]``): the submitted
    request without its configuration, the ``(fit, synth)`` window pairs, and ``extra``
    (``use_linelist``, the line-mask name, a fit's :meth:`FitTab.collect` values)."""
    request = json.loads(request.to_json())
    del request["config"]
    return {"request": request,
            "windows": [[list(fit), list(synth)] for fit, synth in windows], **extra}


class SynthModal:
    def __init__(self, ctx, add, grids, atoms):
        self.ctx = ctx
        self.canvas = ctx.canvas
        self._add = add          # (spec, name, meta=None, conv=None) -> add to the canvas
        ctx.task_importers["synthesis.fit"] = self._import_fit  # a fit another document submitted
        self.grids = grids       # the app-global registry (Preferences), read-only here
        self.atoms = atoms       # the shared NLTE model-atom selection, shown in this dialog
        self.job = None          # ("range", None) | ("masks", LineMaskSet) once opened
        self.source = None       # the name of the spectrum a rerun was opened from
        self.mode = "synth"      # "synth" | "fit"
        self._pending = False    # a submission is validated and closing/closed but the
                                  # actual (slow) handoff to the daemon hasn't run yet

        self.dialog = ctx.bridges.dialog("synth-modal", name="synth-modal-bridge")
        self.title_bridge = ctx.bridges.html("synth-modal-title")
        self.context_bridge = ctx.bridges.html("synth-context")
        self.actions_bridge = ctx.bridges.html("synth-modal-actions")
        self.fit_mode_bridge = ctx.bridges.cls("#synth-modal", "fit-mode")

        # the wavelength window of a range job (Å); seeded from the range that opened the
        # modal, blank from the landing page, and required before a run starts
        self.window_form = Form("synthwindow", [
            Field("xmin", "λ min [Å]", kind="number", value=None, min=0, inline=True,
                  required=True),
            Field("xmax", "λ max [Å]", kind="number", value=None, min=0, inline=True,
                  required=True),
        ])

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

        # the Fit section: target + what varies; the guess is the setup itself
        self.fit = FitTab()
        self.fit.set_catalog(load_prefs())
        self.fit.set_relative(bool(self.setup.base_form.value("relative_abundances")))
        self.fit.set_feh(self.setup.run_form.value("feh"))
        self.fit.targets = [(o.uid, o.name) for o in self.canvas.spectra]
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
        return [*self.setup.models, *self.window_form.models, *self.fit.models]

    # --- opening ------------------------------------------------------------------------
    def open_range(self, xmin=None, xmax=None, mode="synth"):
        """Open on a free window (Å); ``None`` leaves the fields blank for the user."""
        self.job = ("range", None)
        self.source = None
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
        self.source = None
        self.setup.element = mset.element
        self._open(mode)

    def open_rerun(self, obj):
        """Open on the request a synthesis or fit spectrum ``obj`` recorded
        (``meta["rerun"]``, see :func:`rerun_record`): its windows, atmosphere, star,
        abundances, RT mode and model atoms, and for a fit its target and axes, as a draft
        of the setup. The configuration is the current Preferences. Submitting adds a new
        spectrum; ``obj`` is left as it is."""
        record = ((obj.meta or {}).get("rerun") or {}) if obj is not None else {}
        if not record:
            self.notify("This spectrum does not record the request it was computed from.")
            return
        try:
            request = SynthRequest.from_dict({**record["request"],
                                              "config": synth_config(load_prefs()).to_dict()})
            windows = [(tuple(fit), tuple(synth)) for fit, synth in record["windows"]]
        except (KeyError, TypeError, ValueError, OSError) as exc:
            self.notify(f"Cannot read the request of {obj.name}: {exc}")
            return
        if record.get("mask") or len(windows) > 1:
            # the recorded windows, as a set of their own: the source set may have changed
            mset = LineMaskSet(label=record.get("mask") or obj.name, masks=[
                LineMask(lmin=f[0], lmax=f[1], smin=s[0], smax=s[1]) for f, s in windows])
            self.job = ("masks", mset)
        else:
            self.job = ("range", None)
            self.window_form.set_values(xmin=windows[0][0][0], xmax=windows[0][0][1])
        self.setup.element = ""
        self.source = obj.name

        def prefill():
            reason = self.setup.load_request(request, record.get("use_linelist", True))
            self.atoms.selected = list(naming.atoms(request))
            self.atoms.render()
            self.fit.set_relative(bool(self.setup.base_form.value("relative_abundances")))
            self.fit.set_feh(self.setup.run_form.value("feh"))
            if record.get("fit"):
                self.fit.load(record["fit"])
            if reason:
                self.notify(reason)
        self._open("fit" if record.get("fit") else "synth", prefill)

    def _open(self, mode, prefill=None):
        self.draft.begin()
        if prefill is not None:
            prefill()
        self._pending = False   # a fresh job; any earlier submission has long since handed off
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
        source = (f'<div class="hint">Rerun from <b>{html.escape(self.source)}</b>, as '
                  f'recorded; submitting adds a new spectrum.</div>' if self.source else "")
        return source + self._window_html()

    def _window_html(self):
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
                '<div class="hint">Last-used atmosphere setup. Syntheses and fits keep running after the tab closes.</div>'
                f'<div id="{self.window_form.slot_id}" class="synth-window">'
                f'{self.window_form.html()}</div>')

    def actions_html(self):
        if self.mode == "fit":
            return '<button class="primary-btn synth-submit">Fit</button>'
        return '<button class="primary-btn synth-submit">Synthesise</button>'

    def submit(self):
        """The footer button: what the modal was opened for.

        Ignored while an earlier click's submission is still being handed off (see
        :meth:`_begin_submit`) — the button-side click handler already disables itself,
        this is the belt-and-braces guard against a second ``synth-submit`` action that
        reaches the server before that round trip closes it."""
        if self._pending:
            return
        if self.mode == "fit":
            self.run_fit()
        else:
            self.synthesize()

    def _begin_submit(self):
        """Validation has passed: close the dialog and mark the submission in flight
        *before* doing the slow part (JSON-encoding the request and handing it to the
        daemon over a socket). Returning from the action callback now, instead of after
        that round trip, lets Bokeh flush the dialog-close patch to the browser at once —
        closing the modal is the user's feedback that the click was received, not a side
        effect of the submission having finished."""
        self._pending = True
        self.draft.commit()
        self.fit.commit()
        pulse(self.dialog, "close")

    def _end_submit(self, fn):
        """Run ``fn`` (the actual ``ctx.submit_operation`` call) on the next document
        tick, after the close patch above has already gone out, and clear the in-flight
        guard once it returns."""
        def run():
            try:
                fn()
            finally:
                self._pending = False
        self.ctx.schedule(run)

    def _reject(self, message):
        """Refuse the current submit attempt: report why, and reset the busy submit
        button — the dialog stays open, so unlike a successful submit (:meth:`_begin_submit`)
        nothing re-renders the footer or closes the dialog to do that on its own."""
        self.notify(message)
        pulse(self.dialog, "reject")

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
    def _request_rows(cls, request, grid):
        """:meth:`_atmosphere_rows` for a request's atmosphere; ``grid`` is the setup's
        active grid record, which names the grid (or the grid a picked file came from)."""
        atmosphere = request.atmosphere
        if isinstance(atmosphere, AtmosphereGrid):
            return cls._atmosphere_rows({**(grid or {}), "format": atmosphere.format})
        return cls._atmosphere_rows({**(grid or {}), "format": atmosphere.format},
                                    atmosphere.path,
                                    atmosphere.resolution if isinstance(atmosphere, Snapshot)
                                    else None)

    @staticmethod
    def _request_terms(request):
        """``(mode, abundances, relative, atoms)`` of a request, as the Info tab states them."""
        abundances = request.abundances
        text = " ".join(f"{el}={value:g}" for el, value in (abundances or {}).items())
        relative = abundances.relative if abundances is not None else True
        return ("NLTE" if request.is_nlte else "LTE"), text, relative, naming.atoms(request)

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
        self.fit.set_targets([(o.uid, o.name) for o in self.canvas.spectra])

    def run_fit(self):
        """Fit the chosen target as a durable ``synthesis.fit`` workflow.

        The request is the setup's, with the job's windows; the operation refuses what
        fitting does not support (a snapshot, saved departures, stellar axes with a fixed
        file) before the job is admitted."""
        request = self.setup.context()
        if request is None:
            self._reject(self.setup.reason or "Cannot fit — check Preferences → Synthesis.")
            return
        windows = self._job_windows()
        if not windows:
            self._reject("Set the wavelength window (λ min/max) first.")
            return
        c = self.fit.collect()
        obj = self.canvas.spectra.find(c["target"]) if c["target"] else None
        if obj is None:
            self._reject("Select a target spectrum to fit.")
            return
        try:
            request = request.update(windows=[synth for _fit, synth in windows])
            check_fit(request, fit_stellar=c["fit_stellar"], abund_fits=c["abund_fits"],
                      fit_post=c["fit_post"])
        except ValueError as exc:
            self._reject(f"Cannot fit: {exc}")
            return
        grid = self.setup.active_grid()
        star = request.star
        guesses = dict(star)
        mode, abundances, relative, atoms = self._request_terms(request)
        xmin = min(lo for lo, _hi in request.windows)
        xmax = max(hi for _lo, hi in request.windows)
        kind, mset = self.job or ("range", None)
        name = self.ctx.workspace.unique_name(naming.fit(
            obj.name, request, grid, window=mset.name if kind == "masks" else None,
            fitted=c["abund_fits"]))
        rows = self._request_rows(request, grid)
        preparation = self.canvas.preparation(obj)
        rerun = rerun_record(request, windows, mask=mset.name if kind == "masks" else None,
                             fit=c)

        def fit_meta(stellar, terms):
            """The Info payload of the fit: the synthesis, plus the target and the
            preparation it was measured under."""
            meta = self._synth_meta(stellar, mode=mode, abundances=terms, relative=relative,
                                    atoms=atoms, xmin=xmin, xmax=xmax, atmosphere=rows)
            meta["params"].append(["Fit target", " · ".join((
                obj.name, preparation["normalisation"],
                "filtered" if preparation["filtered"] else "unfiltered"))])
            meta["preparation"] = preparation
            meta["rerun"] = rerun
            return meta

        def on_done(result):
            best = result.best
            self.fit.populate(best)
            # the best-fit stellar parameters become the setup, so the next synthesis (or
            # fit) starts from them; the broadening (post) is applied as a live convolution,
            # shown in the Edit tab.
            stellar = {k: float(best.get(k, guesses[k])) for k in guesses}
            self.adopt_stellar(stellar)
            terms = synth.fitted_abundance_terms(
                best, c["abund_fits"], base=abundances, relative=relative,
                feh=float(stellar["feh"]), solar=synth.abund_table(request.config.to_prefs()))
            self._add(result.spectrum, name, meta=fit_meta(stellar, terms), conv=result.post)

        self._begin_submit()
        # the fit only reads the target spectrum inside the synth span, so the record
        # carries that part as measured (filtered when filtering is active), with the
        # guesses as the annotations a document that imports the record later shows, and
        # the target's preparation
        meta = fit_meta(guesses, abundances)
        kwargs = {"target": self.canvas.target_spectrum(obj.uid).crop(xmin=xmin, xmax=xmax),
                  "fit_windows": [list(fit) for fit, _synth in windows],
                  "fit_stellar": list(c["fit_stellar"]), "abund_fits": dict(c["abund_fits"]),
                  "fit_post": dict(c["fit_post"])}
        self._end_submit(lambda: self.ctx.submit_operation(
            name, "synthesis.fit", request, on_done, kwargs=kwargs, meta=meta))

    def _import_fit(self, result, handle):
        """A fit record another document submitted: its best-fit spectrum is added under
        the recorded name and annotations, with the broadening as a live convolution."""
        self.fit.populate(result.best)
        self._add(result.spectrum, handle.label, meta=handle.meta, conv=result.post)

    def _job_windows(self):
        """The job's windows as ``(fit, synth)`` pairs: a mask set's included masks with
        their own fit and synth windows (a mask without synth bounds spans its fit window),
        or the free range as one window with coinciding bounds. Empty when the range fields
        do not hold a valid interval."""
        kind, mset = self.job or ("range", None)
        if kind == "masks":
            windows = []
            for m in mset.masks:
                if not m.include:
                    continue
                synth_window = ((m.smin, m.smax) if m.smin == m.smin and m.smax == m.smax
                                else (m.lmin, m.lmax))
                windows.append(((m.lmin, m.lmax), synth_window))
            return tuple(windows)
        xmin, xmax = self._window_fields()
        return () if xmin is None else (((xmin, xmax), (xmin, xmax)),)

    def adopt_stellar(self, stellar):
        """Write ``{teff, logg, feh, vmic}`` into the setup (form + store)."""
        self.setup.set_stellar(stellar)
        for k, v in stellar.items():
            store.update_pref("synth_run", k, v)
        self.fit.set_feh(stellar["feh"])

    # --- keeping the Fit section in step with the setup --------------------
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

    # --- launch -----------------------------------------------------------
    def synthesize(self):
        """Submit the job as one durable synthesis workflow, whatever its atmosphere."""
        request = self.setup.context()
        if request is None:
            self._reject(self.setup.reason
                        or "Cannot synthesize — check Preferences → Synthesis.")
            return
        windows = self._job_windows()
        if not windows:
            self._reject("Set the synthesis wavelength window (λ min/max) first.")
            return
        try:
            request = request.update(windows=[synth for _fit, synth in windows])
        except ValueError as exc:
            self._reject(f"Cannot synthesize: {exc}")
            return
        kind, mset = self.job or ("range", None)
        star = request.star
        mode, abundances, relative, atoms = self._request_terms(request)
        grid = self.setup.active_grid()
        name = self.ctx.workspace.unique_name(naming.synthesis(
            request, grid, window=mset.name if kind == "masks" else None))
        shown_min = min(lo for lo, _hi in request.windows)
        shown_max = max(hi for _lo, hi in request.windows)
        atmosphere = self._request_rows(request, grid)
        if isinstance(request.departures, Departures):
            atmosphere.append(["Saved departures", str(request.departures.folder)])
        meta = self._synth_meta(dict(star), mode=mode, abundances=abundances,
                                relative=relative, atoms=atoms,
                                xmin=shown_min, xmax=shown_max, atmosphere=atmosphere)
        use_linelist = self.setup.use_linelist
        linelist = ("Off (model atom only)" if not use_linelist else
                    str(request.config.linelist or ""))
        meta["params"].append(["Line list", linelist])
        if kind == "masks":
            meta["params"].append(["Line mask", mset.name])
        meta["rerun"] = rerun_record(request, windows, use_linelist=use_linelist,
                                     mask=mset.name if kind == "masks" else None)

        def add_result(result):
            for index, spec in enumerate(result.spectra, start=1):
                suffix = f"·{index}" if len(result.spectra) > 1 else ""
                self._add(spec, name + suffix, meta=meta)

        self._begin_submit()
        self._end_submit(lambda: self.ctx.submit_operation(
            name, "synthesis.synthesize", request, add_result,
            kwargs={"background_lines": use_linelist}, meta=meta))
