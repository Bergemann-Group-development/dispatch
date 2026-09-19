"""The Ionisation-balance tool: fit stellar parameters to a measured line-mask set.

The workflow starts in the LineMasks view, not here. A line-mask set carrying an element is
the input, and its sidebar's **Fit ionisation balance** action opens this tool's modal: pick
the target observation, an interpolated atmosphere grid, an entry point (a reference-star class or
typed-in numbers), and which of (Teff, log g, [Fe/H], ξ) are fitted rather than held. What
comes back is a :class:`measurement.Measurement` — its own object in the left sidebar,
because everything in it is spectrum-dependent — holding the measured lines and the fit that
stands on them.

Two phases, one durable ``sparams.ionbalance`` workflow, continuing the measurement object
the modal creates (``measured=``) so its snapshots stream into it:

1. **Measure + identify.** :func:`tumult.sparams.measure_balance_lines`: every switched-on
   window's EW, FWHM and EW error on the target, with the dominant transition in each window
   taken from HOP line strengths computed on the MARCS model nearest the entry point
   (:func:`tumult.linestrength.band_on_atmosphere`) rather than a raw log gf ranking — which
   stage a window belongs to is a property of the star, and the ionization balance reads it off this.
   A window another species dominates is flagged ``off_element`` and switched off, and a
   window whose dominant transition is the element's third stage or above simply carries no
   diagnostic weight.
2. **Solve.** :func:`tumult.sparams.ionbalance`: from the entry point, invert every line's
   EW to A(X) through an adaptive curve of growth, compute the four balance diagnostics and
   their Jacobian, and step the *free* parameters (:func:`tumult.sparams.solve_coupled`)
   until the diagnostics are all inside tolerance.

Both phases are headless — this module supplies the modal, the streaming table and the
diagnostic panels, and passes a ``progress`` callback so rows, the panels and the iteration
status line land as they are computed. The view opens on **Balance**, where the two slopes
the solve fits are drawn side by side (A(X) against excitation energy and against reduced EW)
with the two offsets stated as insets, so a fit is watched as it converges; **Lines** shows
the target spectrum with every inverted line's synthetic drawn on it. A held parameter's
diagnostic is still computed and shown, as the residual it is.

Refitting from inside the object reuses its lines and overwrites the fit; refitting a
duplicate (:meth:`resultview.ResultView.duplicate_result`) keeps both fits; launching again from
the LineMasks view measures afresh into a new object. That is the whole identity rule:
the object owns the measurement, the fit is replaceable.
"""
from __future__ import annotations

import html
import time
import logging

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, Label, Slope

from tumult import abundance as ab
from tumult import measurement as ms
from tumult import refstars
from tumult import sparams as sp
from ..workspace import naming
from ..workspace import store
from tumult.synthesis import workflows as synth
from .diagnostics import PANEL_AXES, empty_points
from .localcontinuum import LOCAL_CONTINUUM_HEAD, continuum_cell
from .resultview import ResultView, format_value, preparation_note, value_label
from ..ui.style import ACCENT, MUTED
from ..ui.controls import Field, Form
from ..spectrum.ew import FIT_PROFILES, profile_options
from ..ui.htmlbridge import pulse
from ..preferences.synthprefs import atmosphere_grid, load_prefs, missing_prefs, synth_config
from tumult.synthesis import Star, SynthRequest

ION2_COLOR = "#e05c5c"     # singly-ionised points (the neutral stage uses the shared accent)
# entry-point choice that is not a reference-star class
ENTRY_MANUAL = "manual"

DEFAULTS = {
    "profile": "gauss", "snr": 100.0, "cont_pct": 0.5,
    "entry": refstars.AUTO_KEY,
    "teff": 5750.0, "logg": 4.5, "feh": 0.0, "vmic": 1.0,
    "fit_teff": True, "fit_logg": True, "fit_feh": True, "fit_vmic": True,
    "max_iter": 8.0, "limit": 1.0, "max_rounds": 8.0,
    "tol_exc": sp.DEFAULT_TOLS["s_exc"], "tol_rew": sp.DEFAULT_TOLS["s_rew"],
    "tol_ion": sp.DEFAULT_TOLS["d_ion"], "tol_feh": sp.DEFAULT_TOLS["d_feh"],
}

PARAM_LABELS = {"teff": "Teff", "logg": "logg", "feh": "[Fe/H]", "vmic": "ξ"}
DIAG_LABELS = {"s_exc": "E_low slope", "s_rew": "REW slope",
               "d_ion": "II−I", "d_feh": "metallicity"}
TOL_FIELDS = {"s_exc": "tol_exc", "s_rew": "tol_rew",
              "d_ion": "tol_ion", "d_feh": "tol_feh"}


class BalanceSynthesis:
    """The balance fit's private, interpolation-only synthesis configuration.

    Atmosphere-grid registration and the physical input paths remain application preferences,
    but the selected grid is owned by the balance modal.  Deliberately do not expose an
    individual-atmosphere or pick-file mode: Teff/log g/[Fe/H] cannot be inferred while the
    atmospheric structure is fixed.
    """

    def __init__(self, grids, on_change=None):
        self.grids = grids
        self.on_change = on_change
        self.reason = None
        saved = store.get_section("ionbal")
        ids = [gid for gid, _ in grids.names]
        preferred = saved.get("atmos_grid_id")
        if preferred not in ids:
            preferred = (store.get_section("synthesis").get("ref_grid_id")
                         if store.get_section("synthesis").get("ref_grid_id") in ids
                         else (ids[0] if ids else ""))
        # the curve-of-growth ladders synthesise the element's lines alone unless background
        # lines are included; off by default, as the balance was calibrated that way
        self.form = Form("ionbalsetup", [
            Field("atmos_grid_id", "Atmosphere grid", kind="select", value=preferred,
                  options=grids.names),
            Field("background_lines", "Include background lines", kind="checkbox",
                  value=bool(saved.get("background_lines", False))),
        ], on_change=self._changed)

    @property
    def models(self):
        return self.form.models

    def _changed(self, key, value):
        if self.on_change:
            self.on_change()

    def refresh(self):
        options = list(self.grids.names)
        ids = [gid for gid, _ in options]
        current = str(self.form.value("atmos_grid_id") or "")
        if current not in ids:
            current = ids[0] if ids else ""
        self.form.fields["atmos_grid_id"].options = options or [(
            "", "— no registered atmosphere grid —")]
        self.form.set_values(atmos_grid_id=current)

    def selected_grid(self):
        return self.grids._find(str(self.form.value("atmos_grid_id") or ""))

    def context(self, params):
        """The balance's :class:`~tumult.synthesis.SynthRequest` starting at ``params``
        (``{teff, logg, feh, vmic}``), or ``None`` with :attr:`reason` set."""
        prefs = load_prefs()
        if missing_prefs(prefs):
            self.reason = "Synthesis paths are missing — check Preferences → Synthesis."
            return None
        grid = self.selected_grid()
        if grid is None:
            self.reason = ("Download or register an atmosphere grid in Preferences → Atmosphere grids, then "
                           "select it in this dialog.")
            return None
        if not str(grid.get("folder") or "").strip():
            self.reason = "The selected atmosphere grid has no folder."
            return None
        # A balance is intentionally a reproducible LTE, scaled-solar calculation on an
        # interpolated grid; the element abundance itself is supplied by each
        # curve-of-growth point, and each iteration moves the star.
        reference = self.grids.reference_grid()
        try:
            request = SynthRequest(
                config=synth_config(prefs),
                atmosphere=atmosphere_grid(grid, atmosphere_grid(reference) if reference else None),
                star=Star(**params))
        except (TypeError, ValueError) as exc:
            self.reason = str(exc)
            return None
        self.reason = None
        return request

    @property
    def background_lines(self):
        return bool(self.form.value("background_lines"))


class IonBalanceTool(ResultView):
    """The tool behind the measurement objects: the modal, the solve, and the object view."""

    name = "ionbal"
    prefix = "ionbal"
    kind = "ionbal"
    section = "ionbal"
    steps = 2                        # lines, balance
    empty_table_msg = "No measurement selected"
    paint_overlays = True            # every inverted line is drawn on the target spectrum
    local_continuum = True           # each measured line may carry its own continuum
    y_axis_label = "A(X)"
    final_stage = ms.STAGE_SOLVED
    check_cols = ("include", "continuum")
    extra_sort_keys = {"blend": lambda ln: ln.blend_frac}
    extra_sort_asc = {"blend": True}

    def __init__(self, ctx, *, grids, list_observations, get_observation, master_path):
        self.list_observations = list_observations   # () -> [(uid, name)]
        self.get_observation = get_observation       # (uid) -> Spectrum (normalised) | None
        self.master_path = master_path

        saved = dict(store.get_section(self.section))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # Older sessions offered the Synthesize tab as an entry point.  Preserve the
        # numbers it filled, but make their provenance explicit now that the dependency is
        # gone.
        if vals["entry"] == "setup":
            vals["entry"] = ENTRY_MANUAL

        # --- the modal's three forms: what is measured, where the solve starts, how it runs
        self.form1 = Form("ionbal1", [
            Field("target", "Target spectrum", kind="select", value="", options=[]),
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=profile_options(FIT_PROFILES)),
            Field("snr", "SNR (observation)", kind="number", value=vals["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=vals["cont_pct"], min=0, step=0.1),
        ], on_change=self._on_field, groups=[("Measurement uncertainty", ("snr", "cont_pct"))])
        self.form2 = Form("ionbal2", [
            Field("entry", "Preset", kind="select", value=vals["entry"],
                  options=self._entry_options()),
            Field("teff", "Teff [K]", kind="number", value=vals["teff"], step=50,
                  inline=True),
            Field("logg", "log g", kind="number", value=vals["logg"], step=0.1,
                  inline=True),
            Field("feh", "[Fe/H]", kind="number", value=vals["feh"], step=0.1,
                  inline=True),
            Field("vmic", "Microturbulence ξ [km/s]", kind="number", value=vals["vmic"], step=0.1,
                  inline=True),
            Field("fit_teff", "Fit Teff", kind="checkbox", value=vals["fit_teff"]),
            Field("fit_logg", "Fit log g", kind="checkbox", value=vals["fit_logg"]),
            Field("fit_feh", "Fit [Fe/H]", kind="checkbox", value=vals["fit_feh"]),
            Field("fit_vmic", "Fit ξ", kind="checkbox", value=vals["fit_vmic"]),
        ], on_change=self._on_field)
        self.formadv = Form("ionbaladv", [
            Field("max_iter", "Max iterations", kind="number", value=vals["max_iter"],
                  min=1, step=1),
            Field("limit", "Curve-of-growth half-range [dex]", kind="number",
                  value=vals["limit"], min=0.05, step=0.05),
            Field("max_rounds", "Curve-of-growth rounds", kind="number",
                  value=vals["max_rounds"], min=1, step=1),
            Field("tol_exc", "Tol: excitation slope [dex/eV]", kind="number",
                  value=vals["tol_exc"], min=0, step=0.001),
            Field("tol_rew", "Tol: REW slope [dex/dex]", kind="number",
                  value=vals["tol_rew"], min=0, step=0.005),
            Field("tol_ion", "Tol: ionisation ΔA [dex]", kind="number",
                  value=vals["tol_ion"], min=0, step=0.005),
            Field("tol_feh", "Tol: metallicity consistency [dex]", kind="number",
                  value=vals["tol_feh"], min=0, step=0.005),
        ], on_change=self._on_field)
        setup = BalanceSynthesis(grids, on_change=self._push_modal_context)
        super().__init__(ctx, setup=setup)
        self.own_panel("ionbal-panel-content")
        self._modal = ctx.bridges.dialog("ionbal-modal", name="ionbal-modal-bridge")
        self._context_bridge = ctx.bridges.html_all(".ionbal-context")
        self._pending_set = None     # the mask set the open modal will measure (None = refit)
        self._refit_uid = None       # the measurement the open modal will refit

        # The solve reads four diagnostics off the same lines at once, so the Balance view
        # shows the two slopes as panels — A(X) against excitation energy and against reduced
        # EW — and states the two offsets, which are levels rather than trends, as an inset on
        # the panel whose parameter they belong with: the II−I offset (log g) beside the
        # excitation slope, the metallicity offset beside the REW slope.
        # Neutral-stage points ride the shared source, the ionised stage a second one in red;
        # both panels get whiskers, and each carries the regression over its own *neutral*
        # points.
        self.build_panels()
        self.panels.add_whiskers(self.panels.src, ACCENT)
        self.src2 = ColumnDataSource(empty_points())
        self.src2.name = f"{self.prefix}_points2"
        self.panels.add_series(self.src2, ION2_COLOR)
        self.panels.add_whiskers(self.src2, ION2_COLOR)
        self.trends, self.insets = [], []
        for fig in self.panels.figs:
            trend = Slope(gradient=0, y_intercept=0, line_color=ACCENT,
                          line_dash="dotted", line_width=1.4)
            trend.visible = False
            inset = Label(x=8, y=8, x_units="screen", y_units="screen", text="",
                          text_color=MUTED, text_font_size="11px", background_fill_color=None)
            fig.add_layout(trend)
            fig.add_layout(inset)
            self.trends.append(trend)
            self.insets.append(inset)

    # --- Tool protocol --------------------------------------------------------------------
    @property
    def models(self):
        # All modal controls are private to this tool; only the registered-grid data behind
        # the atmosphere selector comes from application preferences.
        return [*self.form1.models, *self.form2.models, *self.formadv.models,
                *self.setup.models, self.actions]

    def activate(self):
        """Enter the result view; the fit dialog carries its own synthesis context."""
        self.active = True
        self.push_panel()
        self.push_mode("on")
        self.set_step(self.step)
        if self.selected() is not None and self.is_focused():
            self.push_view("on")
            self._render_decor(self.jumped)

    def refresh(self):
        """Keep the modal's target selector in step with the loaded observations."""
        self._sync_targets()
        self.setup.refresh()

    def update(self):
        self.refresh()
        if self.active:
            self.push_panel()
            self._render_decor(self.jumped)   # the highlights follow a display-unit flip

    # --- the modal ------------------------------------------------------------------------
    @staticmethod
    def _entry_options():
        return ([(ENTRY_MANUAL, "Type them in")]
                + refstars.options("Auto (nearest class to the target)"))

    def refresh_refstars(self):
        """The reference-star catalogue changed: relist the entry points; an entry that
        left the list becomes 'Type them in', keeping the numbers it filled."""
        field = self.form2.fields["entry"]
        field.options = self._entry_options()
        if str(field.value) not in {k for k, _ in field.options}:
            field.value = ENTRY_MANUAL
        self.form2.render()

    def _sync_targets(self):
        """Repaint the target dropdown from the canvas (observations come and go)."""
        obs = self.list_observations()
        options = [(str(uid), name) for uid, name in obs]
        field = self.form1.fields["target"]
        if not any(uid == str(field.value) for uid, _ in options):
            field.value = options[0][0] if options else ""
        field.options = options or [("", "— no spectrum loaded —")]
        self.form1.render()

    def open_for(self, mset):
        """Open the modal to measure and fit ``mset`` (the LineMasks action)."""
        if mset is None or not any(m.include for m in mset.masks):
            self.notify("This line mask has no switched-on windows.")
            return
        if not mset.element:
            self.notify("Set the element these windows diagnose before fitting a balance.")
            return
        self._pending_set, self._refit_uid = mset, None
        self._open_modal()

    update_label = "Refit…"

    def update_result(self):
        self.open_refit()

    def open_refit(self, uid=None):
        """Open the modal to refit an existing measurement in place."""
        res = self.find(uid or self.results.selected)
        if res is None:
            return
        self._pending_set, self._refit_uid = None, res.uid
        # the entry point that produced the standing fit is the pre-fill for the next one
        if res.params0:
            self.form2.set_values(entry=res.entry or ENTRY_MANUAL,
                                  **{k: res.params0[k] for k in sp.PARAM_KEYS
                                     if k in res.params0})
        if res.free:
            self.form2.set_values(**{f"fit_{k}": bool(v) for k, v in res.free.items()})
        self._open_modal()

    def _open_modal(self):
        self.draft.begin()
        self._sync_targets()
        self.setup.refresh()
        # a stored entry point is a *rule* ("the nearest class"), so resolve it against the
        # target now rather than showing whatever numbers it produced last time
        _, spec, _ = self._target()
        params = self._entry_params(str(self.form2.value("entry")), spec)
        if params is not None:
            self.form2.set_values(**params)
        self.form2.render()
        self.formadv.render()
        self._push_modal_context()
        pulse(self._modal, "open")

    def _push_modal_context(self):
        if hasattr(self, "_context_bridge"):
            self._context_bridge.value = (self.context_html() + self.setup_html()
                + self.ctx.canvas.preparation_html(self.form1.value("target")))

    def setup_html(self):
        """Compact statement of the modal-owned, interpolation-only synthesis context."""
        grid = self.setup.selected_grid()
        name = (grid or {}).get("name") or "no atmosphere grid"
        lines = ("with background lines" if self.setup.form.value("background_lines")
                 else "element lines only")
        return (f'<div class="hint">Synthesis: interpolated {html.escape(name)} atmospheres '
                f'· LTE · scaled-solar composition · {lines}.</div>')

    def context_html(self):
        """The modal's header block: what is about to be measured, or refitted."""
        res = self.find(self._refit_uid) if self._refit_uid else None
        if res is not None:
            n1, n2 = res.counts()
            el = html.escape(res.element or "X")
            return (f'<div class="hint">Refitting <b>{html.escape(res.name)}</b> — reusing '
                    f'its {len(res.lines)} measured lines ({n1} {el} I, {n2} {el} II) on '
                    f'{html.escape(res.ref_name)}. The measurement is kept; only the fit is '
                    f'replaced.</div>')
        mset = self._pending_set
        if mset is None:
            return '<div class="hint">Nothing selected to fit.</div>'
        n = sum(1 for m in mset.masks if m.include)
        return (f'<div class="hint">Measuring <b>{html.escape(mset.name)}</b> — {n} '
                f'switched-on window{"" if n == 1 else "s"} of '
                f'<b>{html.escape(mset.element)}</b> on the target spectrum, then fitting '
                f'(Teff, log g, [Fe/H], ξ) to their excitation and ionisation balance. Both '
                f'ionisation stages have to be present among them for the balance to '
                f'close.</div>')

    # --- entry point ----------------------------------------------------------------------
    def _entry_params(self, key, target_spec):
        """``{teff, logg, feh, vmic}`` for entry-point ``key`` (None to keep what is typed)."""
        if key == ENTRY_MANUAL:
            return None
        if key == refstars.AUTO_KEY:
            params = self.ctx.strengths.params_of(target_spec) if target_spec is not None \
                and self.ctx.strengths is not None else None
            star = refstars.nearest(*params) if params else refstars.get(refstars.DEFAULT_KEY)
        else:
            star = refstars.get(key)
        return star.params if star is not None else None

    def _target(self):
        """``(uid, spectrum, error)`` for the modal's chosen target observation."""
        uid = str(self.form1.value("target") or "")
        obs = dict((str(u_), n) for u_, n in self.list_observations())
        if uid not in obs:
            return None, None, "Load and select a target observation to measure EWs on."
        spec = self.get_observation(uid)
        if spec is None:
            return None, None, "That observation is no longer loaded."
        if spec.yunit != u.dimensionless_unscaled:
            return None, None, ("Normalise the target first — the balance needs "
                                "continuum-normalised flux.")
        return uid, spec, None

    # --- state persistence ------------------------------------------------------------------
    def _on_field(self, key, value):
        if key == "entry":
            _, spec, _ = self._target()
            params = self._entry_params(str(value), spec)
            if params is not None:
                self.form2.set_values(**params)
        elif key in sp.PARAM_KEYS and str(self.form2.value("entry")) != ENTRY_MANUAL:
            # typing over a class's numbers is what "manual" means — say so rather than
            # letting the next repaint silently put the class's values back
            self.form2.set_values(entry=ENTRY_MANUAL)
            self._render_scatter()

    def _extra_action(self, kind, arg):
        if kind == "run":
            self.run()
        elif kind == "refit":
            self.open_refit(arg or None)

    # --- session persistence ------------------------------------------------------------
    # the envelope is ResultView's; :mod:`results` builds a saved balance
    # (:class:`tumult.measurement.BalanceResult`) and clears the links — its reference
    # observation and the mask set it was measured from — that a load does not restore.

    # --- shared-machinery hooks --------------------------------------------------------------
    def _toggle_line(self, token):
        """Include-toggle a line; on a fitted measurement, re-aggregate the diagnostics from
        the stored per-line abundances (no re-synthesis)."""
        res = self.selected()
        if res is None:
            return
        try:
            idx = int(token)
        except ValueError:
            return
        if not (0 <= idx < len(res.lines)):
            return
        self.jobs.invalidate(res.uid)
        res.lines[idx].include = not res.lines[idx].include
        res.reaggregate()
        self.push_results()
        self.push_panel()
        self._render_selected()

    def _result_summary(self, res):
        return res.label()

    def select_result(self, uid):
        # the highlights are painted after the base implementation, which hands the centre view
        # (and with it the sidebar focus) to this result only at the very end
        super().select_result(uid)
        self._render_decor(self.jumped)

    def _render_selected(self):
        self.push_panel()
        super()._render_selected()
        self._render_decor(self.jumped)

    def _jump_decor(self, res, ln):
        # the windows and the inverted syntheses are already on the canvas — jumping hatches
        # this line's measured profile and moves the accent onto it.
        self._set_hatch(res, ln)
        self.jumped = next((i for i, l in enumerate(res.lines) if l is ln), None)
        self._render_decor(self.jumped)

    def _result_masks(self, res):
        """The measurement's windows: each LineResult embeds the mask it was measured over,
        so these *are* the windows the balance was solved on.

        Until it has its lines (the view opens the moment the fit is submitted, and measuring
        runs on the worker) they are taken straight from the source mask set, whose own
        switched-off masks are exactly the ones this run skipped."""
        if res.lines:
            return super()._result_masks(res)
        mset = self._source_set(res)
        return list(mset.masks) if mset is not None else []

    def _result_spec(self, res):
        return self._spec_of(res)

    def _clear_line_decor(self):
        """Drop the per-line overlay/hatch and the accent; the highlights stay while a result
        is in view (``_render_decor`` clears them once it is not)."""
        self.jumped = None
        super()._clear_line_decor()
        self._render_decor()

    def _spec_of(self, res):
        """The reference observation, by identity alone (``None`` once it is gone)."""
        return self.get_observation(res.ref_uid) if res.ref_uid is not None else None

    def _source_set(self, res):
        """The live mask set this measurement came from, or None once it is gone.

        ``maskset_uid`` is stored as a string (it survives a session snapshot as JSON) while
        the store keys its sets by int, so the match is made on the string form."""
        masks = getattr(self.ctx, "masks", None)
        if masks is None or not res.maskset_uid:
            return None
        return next((m for m in masks.sets if str(m.uid) == str(res.maskset_uid)), None)

    def on_filter_changed(self):
        """Re-measure the selected measurement's EWs when the canvas Filtered switch toggles,
        with the SNR and continuum uncertainty the measurement was made with. Only before a
        fit stands on the lines: a solved result's abundances were inverted against the EWs
        measured then, so it is marked out of date instead and a refit measures them anew
        (:meth:`sync_target`)."""
        if not self.active:
            return
        res = self.selected()
        if res is None:
            return
        if res.stage != ms.STAGE_MEASURED:
            self.push_results()
            self._render_selected()
            return
        spec = self._spec_of(res)
        if spec is None or spec.yunit != u.dimensionless_unscaled:
            return
        outcome = ab.refresh_ews(spec, res.lines, profile=res.profile, snr=res.snr,
                                 cont_pct=res.cont_pct)
        self.report_refresh(res, outcome)

    # --- the run ------------------------------------------------------------------------------
    def run(self):
        """The modal's action: measure (or reuse) the lines, then solve."""
        refit = self.find(self._refit_uid) if self._refit_uid else None
        mset = self._pending_set
        if refit is None and mset is None:
            self.notify("Nothing to fit — start from a line mask.")
            pulse(self._modal, "reject")
            return
        ref_uid, spec, err = self._target() if refit is None else \
            (None, self._spec_of(refit), None)
        if refit is not None and spec is None:
            err = "The target observation for this measurement is no longer loaded."
        if err:
            self.notify(err)
            pulse(self._modal, "reject")
            return
        v2 = self.form2.values()
        free = {k: bool(v2[f"fit_{k}"]) for k in sp.PARAM_KEYS}
        entry = str(v2["entry"])
        params0 = {k: float(v2[k]) for k in sp.PARAM_KEYS}
        ctx = self.setup.context(params0)
        if ctx is None:
            self.notify(self.setup.reason or "Configure the synthesis setup first.")
            pulse(self._modal, "reject")
            return
        va = self.formadv.values()
        tols = {k: float(va[f]) for k, f in TOL_FIELDS.items()}

        if refit is None:
            v1 = self.form1.values()
            res = ms.BalanceResult(
                element=mset.element, ref_uid=ref_uid, maskset_uid=str(mset.uid), maskset_name=mset.name,
                maskset_key=ms.mask_key(mset), ref_name=dict(self.list_observations())
                .get(ref_uid, str(ref_uid)),
                profile=str(v1["profile"]), snr=float(v1["snr"]) or None,
                cont_pct=float(v1["cont_pct"]))
            res.preparation = self.ctx.canvas.preparation(self.ctx.canvas.spectra.find(ref_uid))
            res.preparation["line_list"] = str(self.master_path)
            self._add_result(res, naming.analysis(
                res.element, "balance", res.ref_name, mset.name, ctx,
                self.setup.selected_grid()))
        else:
            res = refit
            self.sync_target(res, spec)
        res.params0, res.entry, res.free, res.tols = dict(params0), entry, free, tols
        res.stage = ms.STAGE_MEASURED
        res.history, res.diags, res.sigmas, res.converged = [], {}, {}, False
        # the progress strip shows from the submission on; the record's snapshots take over
        res.progress = {"stage": "Measuring windows" if refit is None else "Screening lines",
                        "done": 0, "total": 1 + sp.EXPECTED_ITERATIONS, "started": time.time()}

        self.draft.commit()
        pulse(self._modal, "close")
        self.push_results()
        self.set_step(1)             # land on the balance panels: the fit is what to watch
        self.select_result(res.uid)
        self.notify(f"{res.element} balance: "
                    + ("solving…" if refit is not None else "measuring the windows…"))
        kwargs = {"target": spec, "measured": res, "free": free,
                  "tols": tols, "max_iter": max(1, int(float(va["max_iter"]))),
                  "limit": float(va["limit"]),
                  "max_rounds": max(1, int(float(va["max_rounds"]))),
                  "background_lines": self.setup.background_lines}
        if refit is None:
            kwargs["line_mask"] = mset
        # the close patch above only reaches the browser once this callback returns, so the
        # actual submission - a blocking round trip to the workflow daemon - is deferred one
        # document tick (as synthmodal.SynthModal._begin_submit/_end_submit does), rather than
        # sitting in front of it and leaving the modal looking inert until it answers
        self.ctx.schedule(lambda: self.submit(res, res.name, self.task_kind, ctx,
                                              self._on_done, kwargs=kwargs))

    task_kind = "sparams.ionbalance"
    import_task_result = lambda self, result, handle: self._on_done(
        self._adopt(result, handle.label))

    def _on_done(self, result):
        self.ctx.emit("masks-changed")
        result.stage = ms.STAGE_SOLVED
        self.push_results()
        self.push_panel()
        if self.selected() is result:
            self._render_selected()
        p, s = result.params, result.sigmas
        if result.diagnose_only:
            self.notify(f"{result.element} balance diagnostics at the given parameters — "
                        f"see the residuals under the table.", level="info")
            return
        if not any(np.isfinite((result.diags or {}).get(k, float("nan")))
                   for k in sp.DIAG_KEYS):
            self.notify(f"No {result.element} line could be inverted at these parameters — "
                        f"the flags say why (out_of_range means the entry point is too far "
                        f"from the star for the curve of growth to reach the measured EWs).")
            return
        state = result.status_label
        parts = " ".join(
            f"{PARAM_LABELS[k]}={p[k]:.2f}" if k != "teff" else f"Teff={p[k]:.0f}"
            for k in sp.PARAM_KEYS)
        self.notify(f"{result.element} balance {state} after {len(result.history)} it: {parts}",
                    level="success" if result.converged else "warn")

    # --- rendering --------------------------------------------------------------------------
    # the offset each panel states as its inset, beside the slope it is fitted with
    PANEL_INSETS = {"evi": "d_ion", "rew": "d_feh"}

    def _render_scatter(self):
        """Repaint both diagnostic panels from the lines as they stand.

        Called on every streamed row as well as at the end, so the panels are the live view
        of the solve: the points move as each iteration re-inverts the EWs, the neutral-stage
        regression follows them, and the insets restate the two offsets.
        """
        res = self.selected()
        inc = [(i, ln) for i, ln in enumerate(res.lines)
               if ln.include and np.isfinite(ln.a)] if res else []
        if not inc:
            self.panels.src.data = empty_points()
            self.src2.data = empty_points()
            self.panels.set_mean(None)
            for trend, inset in zip(self.trends, self.insets):
                trend.visible = False
                inset.text = ""
            return
        self.panels.set_y_label(value_label(res))

        st1 = [(i, ln) for i, ln in inc if sp.is_ion_stage(ln.ion, 1)]
        st2 = [(i, ln) for i, ln in inc if sp.is_ion_stage(ln.ion, 2)]
        d1, d2 = self.panels.pack_points(st1), self.panels.pack_points(st2)
        self.panels.src.data = d1
        self.src2.data = d2
        self.panels.crop_panels((d1, d2))

        # each panel's own weighted regression over the neutral-stage points — the slope the
        # solver reads its temperature (E_low) and microturbulence (REW) from
        sig1 = [ln.sigma_a for _i, ln in st1]
        for key, trend in zip(PANEL_AXES, self.trends):
            slope, _, intercept = sp.weighted_slope(d1[f"x_{key}"], d1["a"], sig1)
            if np.isfinite(slope):
                trend.gradient, trend.y_intercept = slope, intercept
                trend.visible = True
            else:
                trend.visible = False
        for key, inset in zip(PANEL_AXES, self.insets):
            inset.text = self._inset_text(res, self.PANEL_INSETS[key])

        diags = res.diags or {}
        self.panels.set_mean(diags.get("a_mean", float("nan")),
                             diags.get("d_feh_err", float("nan")))

    def _inset_text(self, res, key):
        """A panel's inset: the offset diagnostic it states, with its error and pass mark."""
        diags = res.diags or {}
        v, err = diags.get(key, float("nan")), diags.get(f"{key}_err", float("nan"))
        if not np.isfinite(v):
            return ""
        tol = (res.tols or sp.DEFAULT_TOLS).get(key, sp.DEFAULT_TOLS[key])
        param = next(k for k, dk in sp.DIAG_FOR.items() if dk == key)
        mark = "(held)" if not res.free.get(param, True) else \
            ("\u2713" if abs(v) <= tol else "\u2717")
        out = f"{DIAG_LABELS[key]} {v:+.3f}"
        if np.isfinite(err):
            out += f" \u00b1 {err:.3f}"
        return f"{out} {mark}"

    def results_html(self):
        if not self.results:
            return ""
        rows = self._results_rows("B")
        return ('<div class="section-header result-section-head">Ionisation balance</div>'
                + "".join(rows))

    def panel_html(self):
        """The right sidebar while a measurement owns the centre view."""
        res = self.selected()
        if res is None:
            return '<div class="hint">Select a measurement object.</div>'
        el = html.escape(res.element or "X")
        n1, n2 = res.counts()
        mset = self._source_set(res)
        stale = ('<span class="lm-dirty" title="The source line mask has been edited since '
                 'this was measured">Source edited</span>' if res.stale_against(mset) else "")
        src = html.escape(res.maskset_name or "—")
        rows = [f'<div class="lm-title"><span class="section-header">{el} balance</span>'
                f'{stale}</div>',
                self.stale_banner(res),
                f'<div class="hint">{len(res.lines)} windows from {src}, measured on '
                f'{html.escape(res.ref_name)} · {n1} {el} I, {n2} {el} II included.</div>',
                preparation_note(res)]
        if res.lines and not res.has_both_stages():
            rows.append('<div class="hint hint-warn">Only one ionisation stage among the '
                        'included windows — the II−I offset, and with it log g, cannot be '
                        'measured from this set.</div>')
        if res.ident_note:
            rows.append(f'<div class="hint">Stages identified by '
                        f'{html.escape(res.ident_note)}.</div>')
        if res.params0:
            entry = self._entry_label(res.entry)
            held = [PARAM_LABELS[k] for k in sp.PARAM_KEYS if not res.free.get(k, True)]
            rows.append(f'<div class="hint">Started from {html.escape(entry)}'
                        + (f' · held {html.escape(", ".join(held))}' if held else "")
                        + '.</div>')
        rows.append(self._params_html(res))
        rows.append('<button class="text-btn ionbal-refit">Refit…</button>'
                    '<div class="hint">Refitting keeps these measured lines and their '
                    'curation (measuring their EWs again if the target changed), and '
                    'replaces the fit standing on them; refit a copy (Edit → Duplicate) to '
                    'keep this fit. To measure again — a different star, or edited windows '
                    '— start from the LineMasks view.</div>')
        return "".join(rows)

    @staticmethod
    def _entry_label(key):
        if key == "setup":                 # legacy saved measurement
            return "the former synthesis setup"
        if key == ENTRY_MANUAL:
            return "typed-in parameters"
        if key == refstars.AUTO_KEY:
            return "the target's nearest class"
        return refstars.get(key).label

    def _params_html(self, res):
        """The solved (or held) parameters block."""
        if res.stage != ms.STAGE_SOLVED or not res.params:
            return '<div class="hint">Not fitted yet.</div>'
        stale = bool(self.stale_reasons(res))
        cells = []
        for k in sp.PARAM_KEYS:
            digits = 0 if k == "teff" else 2
            val = format_value(res.params.get(k), res.sigmas.get(k), digits=digits,
                            plus=(k == "feh"))
            held = "" if res.free.get(k, True) else ' <span class="lm-flag">held</span>'
            cells.append(f'<div class="field-group field-inline"><label>{PARAM_LABELS[k]}'
                         f'</label><span class="field-static">{val}{held}</span></div>')
        state = res.status_label
        head = ('<div class="section-header">Parameters <span class="stale-tag">'
                'Previous solution</span></div>' if stale else
                '<div class="section-header">Parameters</div>')
        return (f'<div class="params-block{" stale-block" if stale else ""}">{head}'
                f'{"".join(cells)}'
                f'<div class="hint">{state} after {len(res.history)} iteration'
                f'{"" if len(res.history) == 1 else "s"}.</div></div>')

    def _diag_cell(self, res, key, digits=3):
        d = res.diags or {}
        v = d.get(key, float("nan"))
        if not np.isfinite(v):
            return "—"
        tol = (res.tols or sp.DEFAULT_TOLS).get(key, sp.DEFAULT_TOLS[key])
        param = next(k for k, dk in sp.DIAG_FOR.items() if dk == key)
        if not res.free.get(param, True):
            # a held parameter's diagnostic is a residual, not a target: no pass/fail mark
            return f"{v:+.{digits}f} (held)"
        return f"{v:+.{digits}f} {'✓' if abs(v) <= tol else '✗'}"

    def _status_html(self, res):
        """The iteration/convergence block under the table (streams during the solve)."""
        if not res.history:
            return ""
        last = res.history[-1]
        p, d = res.params or last["params"], res.diags or last["diags"]
        el = html.escape(res.element or "X")
        head = (f'iter {len(res.history)} · Teff {p["teff"]:.0f} · logg {p["logg"]:.2f} · '
                f'[Fe/H] {p["feh"]:+.2f} · ξ {p["vmic"]:.2f}')
        cells = " · ".join(
            f'{DIAG_LABELS[k]} {self._diag_cell(res, k)}' for k in sp.DIAG_KEYS)
        cells += f' (n = {d.get("n_1", 0)} {el} I + {d.get("n_2", 0)} {el} II)'
        stale = False
        if res.stage == ms.STAGE_SOLVED:
            s = res.sigmas or {}
            state = res.status_label
            stale = bool(self.stale_reasons(res))
            head = (('<span class="stale-tag">Previous solution</span> ' if stale else "")
                    + f'<b>{state}</b> after {len(res.history)} it · '
                    f'Teff {format_value(res.params.get("teff"), s.get("teff"), digits=0)} K · '
                    f'logg {format_value(res.params.get("logg"), s.get("logg"), digits=2)} · '
                    f'[Fe/H] {format_value(res.params.get("feh"), s.get("feh"), plus=True, digits=2)} · '
                    f'ξ {format_value(res.params.get("vmic"), s.get("vmic"), digits=2)} km/s')
        return (f'<div class="abund-summary ionbal-status{" stale-block" if stale else ""}">'
                f'{head}'
                f'<span class="abund-budget"> — {cells}</span></div>')

    def table_html(self, res):
        el = html.escape(res.element or "X")
        columns = (
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), (LOCAL_CONTINUUM_HEAD, "continuum"), ("blend", "blend"),
            (f"A({el})", "a"), ("flags", None))
        head = self._table_head(res, columns)
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
            chk = " checked" if ln.include else ""
            cls = " row-selected" if self.jumped == i else ""
            cls += "" if ln.include else " line-excluded"
            ew = format_value(ln.ew, ln.ew_err, digits=1)
            blend = format_value(ln.blend_frac, digits=2)
            a = format_value(ln.a, ln.sigma_a if res.stage == ms.STAGE_SOLVED else None,
                          plus=True)
            fl = ", ".join(ln.flags)
            body.append(
                f'<tr class="ionbal-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="line-include-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'{continuum_cell(ln, i)}<td>{blend}</td><td>{a}</td>'
                f'<td class="abund-flags">{html.escape(fl)}</td></tr>')
        if not res.lines:
            return f'<div class="linelist-empty">No {el} lines measured yet.</div>'
        return (f'<table class="linelist-table abund-table ionbal-lines" data-tool="ionbal">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body) or self._filtered_out_row(columns)}</tbody></table>{self._status_html(res)}')
