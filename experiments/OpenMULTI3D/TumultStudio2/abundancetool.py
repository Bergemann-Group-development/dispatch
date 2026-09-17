"""The Abundance tool: an uncertainty-aware **2-step wizard** for element abundances.

The workflow starts in the LineMasks view. A line-mask set carrying an element is the input,
and its sidebar's **Measure abundances** action opens this tool's modal, which holds both
halves of a run: what is measured (the reference observation, the EW profile, the SNR and
continuum-% behind each EW error) and the synthesis it is measured against (the tool's own
:class:`synthsetup.SynthSetup` — atmosphere, stellar parameters, abundances, RT mode). One
submit therefore does both, as the balance modal does: every switched-on window is measured —
EW, FWHM and EW error, with the dominant transition identified from the line list — and each
line's curve of growth is then synthesised and inverted to the abundance reproducing its EW,
**streaming into the table** as each lands.

What comes back is an :class:`abundance.AbundanceResult` object in the left sidebar and a
split view (spectrum on top, per-line table below) in the centre, with the set's windows drawn
on the spectrum exactly as the LineMasks view draws them. The user curates the lines — clicking
a row jumps the figure to that line, accents its window and overlays its best-fit synthetic.

The right-sidebar panel then carries the two steps that stand on that measurement:

1. **Synthesis**: re-run the per-line synthesis over the curated lines, and reopen the modal
   to change what it runs on. The measurement itself is never redone here: a different
   reference or a different set of windows is a new run, started from the LineMasks view.
   Every line's best-fit synthetic is drawn on the observation as it lands, so the run is
   watched rather than waited on.
2. **Uncertainties**: compute each line's EW and abundance error budget, drawn as error bars
   on the two diagnostic panels that replace the spectrum pane — A(X) against excitation
   energy and against reduced EW, the same pair an ionisation balance is read against, so a
   temperature or microturbulence trend left in the abundances is visible here.

The measurement is :func:`tumult.abundance.measure_lines`, run in this document so the lines
can be curated; the synthesis is the ``abundance.fit`` operation with ``method="profile"``
continuing that measurement (``measured=``) on the setup's request, as a durable workflow
whose progress snapshots stream each line into the table as it lands; the uncertainties are
:func:`tumult.abundance.finalize_result` with the SNR and continuum uncertainty stored on the
result when its windows were measured. What is left in this module is the wizard's own
machinery, shared with the corrections tool through :mod:`abundwizard`, and the rendering.
"""
from __future__ import annotations

from copy import deepcopy

import html
import logging
from pathlib import Path

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource

from tumult import abundance as ab
import store
from tumult.synthesis import workflows as synth
from abundwizard import LOCAL_CONTINUUM_HEAD, PROFILES, AbundanceWizard, preparation_note
from canvas import ACCENT
from controls import Field, Form
from ew import profile_options
from htmlbridge import pulse
import naming

DEFAULTS = {
    "profile": "gauss", "snr": 100.0, "cont_pct": 0.5,
}


class AbundanceTool(AbundanceWizard):
    name = "abundance"
    prefix = "abund"
    kind = "result"
    empty_table_msg = "No abundance result selected"
    paint_overlays = True            # every synthesised line is drawn on the spectrum
    steps = 2                        # synthesis, uncertainties (the lines come from the modal)
    y_axis_label = "A(X)"
    final_stage = "uncert"
    # columns beyond the shared line attributes (see AbundanceWizard.BASE_SORT_KEYS)
    extra_sort_keys = {
        "blend": lambda ln: ln.blend_frac,
        "vshift": lambda ln: ln.vshift,
        "vdop": lambda ln: ln.vdop,
    }
    extra_sort_asc = {"blend": True, "vshift": True, "vdop": True}

    def __init__(self, ctx, *, setup, list_observations, get_observation, master_path):
        self.list_observations = list_observations  # () -> [(uid, name)]
        self.get_observation = get_observation      # (uid) -> Spectrum (normalised) | None
        self.master_path = master_path

        saved = dict(store.get_section("abundance"))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # the modal: how the set handed over by the LineMasks view is measured (independent
        # of stellar parameters). The EW error inputs (SNR / continuum-%) live here too, so
        # each line's EW uncertainty is known as soon as it is measured; a line whose photon
        # EW error exceeds its EW is flagged automatically. The windows themselves are taken
        # exactly as curated — which element they diagnose is the set's business, not this
        # tool's.
        self.form1 = Form("abund1", [
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=profile_options(PROFILES)),
            Field("snr", "SNR (observation)", kind="number", value=vals["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=vals["cont_pct"], min=0, step=0.1),
        ], on_change=self._on_field, groups=[("Measurement uncertainty", ("snr", "cont_pct"))])
        super().__init__(ctx, setup=setup)
        self._ref_bridge = ctx.bridges.html("abund-ref")   # the modal's reference selector
        self._reconfigure = ctx.bridges.cls("#abund-modal", "reconfigure")
        self._modal = ctx.bridges.dialog("abund-modal", name="abund-modal-bridge")
        self._context_bridge = ctx.bridges.html_all(".abund-context")
        self._pending_set = None      # the mask set the open modal will measure and synthesise
        self._settings_uid = None     # ... or the result whose synthesis it reconfigures
        self._actions_bridge = ctx.bridges.html("abund-modal-actions")
        self._note_bridge = ctx.bridges.html("abund-setup-note")   # panel: the run's context
        self._stale_bridge = ctx.bridges.html("abund-stale")       # panel: out-of-date banner
        self.ref_uid = saved.get("ref_uid")

        # the pair of diagnostic panels (excitation and reduced EW), with the per-line
        # abundance error as a bar on each — the same two the balance is read against, so a
        # trend left in either is visible here too
        self._build_figs(ColumnDataSource(self.empty_points()))
        self._add_whiskers(self.src, ACCENT)

    # --- Tool protocol -----------------------------------------------------------------
    def activate(self):
        # the highlights are painted after the base implementation, which hands the centre
        # view (and with it the sidebar focus) to the selected result at the very end
        super().activate()
        self.push_setup_note()
        self._render_decor(self._jumped)

    def refresh(self):
        """Keep the modal's reference selector in step with the loaded observations. Pushed
        whether or not the panel is open: the modal is opened from the LineMasks view."""
        self._ref_bridge.value = self.ref_html()

    def update(self):
        self.refresh()
        if self.active:
            self.push_stale()                  # a continuum edit may have changed the target
            self._render_decor(self._jumped)   # the highlights follow a display-unit flip

    # --- the modal ----------------------------------------------------------------------
    def open_for(self, mset):
        """Open the modal to measure ``mset`` and synthesise it (the LineMasks action)."""
        if mset is None or not any(m.include for m in mset.masks):
            self.notify("This line mask has no switched-on windows.")
            return
        if not mset.element:
            self.notify("Set the element these windows diagnose before measuring them.")
            return
        self._pending_set, self._settings_uid = mset, None
        self.setup.element = mset.element
        self._open_modal()

    def open_settings(self):
        """Open the modal on the selected result to change what its synthesis runs on.

        The measurement half is left in place but does nothing here: those lines were measured
        once, and re-measuring them under different settings would be a different run — which
        is what starting again from the LineMasks view gives.
        """
        res = self.selected()
        if res is None:
            self.notify("Select an abundance object first.")
            return
        self._pending_set, self._settings_uid = None, res.uid
        self.setup.element = res.element
        self._open_modal()

    def _open_modal(self):
        self.draft.begin()
        self.refresh()
        pulse(self._reconfigure, "on" if self._settings_uid else "off")
        self.setup.render()          # the synthesis half lives in the modal
        self._push_modal_context()
        pulse(self._modal, "open")

    def _push_modal_context(self):
        self._context_bridge.value = self.context_html()
        self._actions_bridge.value = self.actions_html()

    def context_html(self):
        """The modal's header block: what is about to be measured, or reconfigured."""
        res = self._find(self._settings_uid) if self._settings_uid else None
        if res is not None:
            return (f'<div class="hint">Settings for <b>{html.escape(res.name)}</b> — its '
                    f'{len(res.lines)} measured lines are kept; synthesising again replaces '
                    f'the abundances standing on them. Duplicate the object first to keep '
                    f'them.</div>')
        mset = self._pending_set
        if mset is None:
            return '<div class="hint">Nothing selected to measure.</div>'
        n = sum(1 for m in mset.masks if m.include)
        return (f'<div class="hint">Measuring <b>{html.escape(mset.name)}</b> — {n} '
                f'switched-on window{"" if n == 1 else "s"} of '
                f'<b>{html.escape(mset.element)}</b> on the reference below, then synthesising '
                f'each line on the atmosphere beside it to invert its EW to '
                f'A({html.escape(mset.element)}) with an error budget. The windows are taken '
                f'as curated; the result is an abundance object of its own.</div>')

    def actions_html(self):
        """The modal's footer: one submit, doing what the modal was opened for."""
        if self._settings_uid:
            return ('<button class="primary-btn abund-run-synth">Synthesise again</button>')
        return '<button class="primary-btn abund-run-lines">Measure &amp; synthesise</button>'

    # --- the object panel's statement of what the synthesis runs on ---------------------
    def push_setup_note(self):
        self._note_bridge.value = self.setup_note_html()

    def setup_note_html(self):
        """Compact statement of the setup the modal owns, so the panel's Run synthesis says
        what it would run on without carrying the fields."""
        grid = self.setup.active_grid() or {}
        name = grid.get("name") or "no atmosphere grid"
        v = self.setup.run_form.values()
        try:
            params = (f"Teff {float(v['teff']):.0f} · log g {float(v['logg']):.2f} · "
                      f"[Fe/H] {float(v['feh']):+.2f} · ξ {float(v['vmic']):.2f} km/s")
        except (TypeError, ValueError):
            params = "stellar parameters unset"
        mode = "NLTE" if self.setup.is_nlte else "LTE"
        file_note = ""
        if self.setup.atmos_file():
            file_note = f" · {html.escape(Path(self.setup.atmos_file()).name)}"
        return (preparation_note(self.selected())
                + f'<div class="hint">Synthesis: {html.escape(name)}{file_note} · {params} · '
                f'{mode}.</div>')

    def on_filter_changed(self):
        """Re-measure the selected result's EWs when the canvas Filtered switch toggles.

        ``get_observation`` already returns the filter-effective reference, so this re-runs the
        shared :func:`tumult.abundance.refresh_ews` over each line's stored fit window with the
        SNR and continuum uncertainty the measurement was made with and repaints the table when
        anything changed. Restricted to the ``"lines"`` stage on purpose: once synthesis has
        run, the abundances were inverted against the EWs measured then, so the result is
        marked out of date instead (:meth:`stale_banner`) and **Synthesise again** measures
        them anew. Windows, flags and include toggles are left untouched."""
        if not self.active:
            return
        res = self.selected()
        if res is None:
            return
        if res.stage != "lines":
            self.push_results()
            self._render_selected()
            return
        spec = self.get_observation(getattr(res, "ref_uid", None))
        if spec is None or spec.yunit != u.dimensionless_unscaled:
            return
        outcome = ab.refresh_ews(spec, res.lines, profile=res.profile, snr=res.snr,
                                 cont_pct=res.cont_pct)
        self.report_refresh(res, outcome)

    # --- session persistence ------------------------------------------------------------
    # the envelope, the factory and the observation relink are AbundanceWizard's: a result
    # names its reference observation by uid alone, cleared when that spectrum is not
    # restored.

    # --- state persistence -------------------------------------------------------------
    def _on_field(self, key, value):
        pass  # The calculation draft is persisted on submission.

    def _extra_action(self, kind, arg):
        if kind == "ref":
            self.ref_uid = arg or None
        elif kind == "run-lines":
            self.run_lines()
        elif kind == "settings":
            self.open_settings()
        elif kind == "run-synth":
            self.run_synth()
        elif kind == "run-uncert":
            self.run_uncert()

    # --- shared-machinery hooks -----------------------------------------------------------
    def _summarise(self, lines):
        return ab.aggregate(lines)

    def _result_summary(self, res):
        s = res.summary
        if np.isfinite(s.get("a_mean", np.nan)):
            return (f"⟨{self.value_label(res)}⟩ = "
                    f"{s['a_mean']:.2f} ± {s['sigma_mean']:.2f}")
        return f"{len(res.lines)} lines measured" if res.stage == "lines" else ""

    # --- window highlights (the set's windows, as the LineMasks view draws them) --------
    def select_result(self, uid):
        # the highlights are painted after the base implementation, which hands the centre
        # view (and with it the sidebar focus) to this result only at the very end
        super().select_result(uid)
        self._render_decor(self._jumped)

    def _render_selected(self):
        super()._render_selected()
        self.push_stale()
        self._render_decor(self._jumped)

    # --- out of date ----------------------------------------------------------------------
    update_label = "Synthesise again"

    def push_stale(self):
        self._stale_bridge.value = self.stale_banner(self.selected())

    def update_result(self):
        """Synthesise the selected result again with the request it was computed from, after
        measuring its EWs again when its target changed."""
        res = self.selected()
        if res is None:
            return
        spec = self.get_observation(res.ref_uid)
        if spec is None:
            self.notify("The reference observation is no longer loaded.")
            return
        self.sync_target(res, spec)
        self.run_synth(request=res.request)

    def _result_spec(self, res):
        """The reference the result's lines were measured on (its points snap the window
        bounds, so window membership is unambiguous)."""
        return self.get_observation(getattr(res, "ref_uid", None))

    def _clear_line_decor(self):
        """Drop the per-line overlay/hatch and the accent; the highlights stay while a result
        is in view (``_render_decor`` clears them once it is not)."""
        self._jumped = None
        super()._clear_line_decor()
        self._render_decor()

    def _jump_decor(self, res, ln):
        # every line's synthetic is already on the canvas — jumping hatches this line's
        # measured EW and moves the accent onto its window
        self.canvas.set_abund_hatch(ln.fit_func, ln.fit_window,
                                        nodes=ln.continuum if ln.renorm else None)
        self._jumped = next((i for i, l in enumerate(res.lines) if l is ln), None)
        self._render_decor(self._jumped)

    # --- reference resolution ----------------------------------------------------------
    def _ref(self):
        """(ref_spec, error) — the selected observation, validated normalised, or an error."""
        obs = self.list_observations()
        ref_uid = self.ref_uid or (obs[0][0] if obs else None)
        spec = self.get_observation(ref_uid) if ref_uid else None
        if spec is None:
            return None, None, "Load and select a reference observation to measure EWs on."
        if spec.yunit != u.dimensionless_unscaled:
            return None, None, ("Normalise the reference first — abundances need "
                                "continuum-normalised flux.")
        return ref_uid, spec, None

    # --- the modal's run: measure the set's windows on the reference --------------------
    def run_lines(self):
        ref_uid, ref_spec, err = self._ref()
        if err:
            self.notify(err)
            return
        mset = self._pending_set
        if mset is None:
            self.notify("Nothing to measure — start from a line mask.")
            return
        if not any(m.include for m in mset.masks):
            self.notify(f"Every window in “{mset.name}” is switched off.")
            return
        # the synthesis runs straight after the measurement, so its setup is checked here —
        # measuring first and failing afterwards would leave a half-run object behind
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason or "Configure the synthesis setup first.")
            return
        v = self.form1.values()
        preparation = self.ctx.canvas.preparation(self.ctx.canvas._find(ref_uid))
        preparation["line_list"] = str(self.master_path)
        self.draft.commit()
        pulse(self._modal, "close")
        self.notify(f"Measuring {mset.element} windows on the reference…")
        mset, ref_spec = deepcopy(mset), deepcopy(ref_spec)
        name = self.ctx.unique_name(naming.analysis(
            mset.element, "abundance", dict(self.list_observations()).get(ref_uid, ""), mset.name,
            request, self.setup.active_grid()))
        loader = lambda _res, token, progress: self._measure_lines(
            ref_uid, ref_spec, mset, str(v["profile"]), float(v["snr"]) or None,
            float(v["cont_pct"]))
        def done(res):
            res.preparation = preparation
            self._on_lines_done(res, name, request, ref_spec)
        self.jobs.submit(f"abund-lines:{name}", loader, done)

    def _measure_lines(self, ref_uid, ref_spec, mset, profile, snr, cont_pct):
        """Measure the set's switched-on windows (:func:`tumult.abundance.measure_lines`)."""
        res = ab.measure_lines(mset, ref_spec, linelist=self.master_path,
                               profile=profile, snr=snr, cont_pct=cont_pct)
        res.ref_uid = ref_uid
        res.ref_name = dict(self.list_observations()).get(ref_uid, "")  # for messages only
        res.snr, res.cont_pct = snr, cont_pct
        return res

    def _on_lines_done(self, result, name, request, observation):
        self._add_result(result, name)
        self.set_step(1)                 # the synthesis is what stands on the measurement now
        self.select_result(result.uid)   # focus it: enters abund-mode, renders the split view
        self.notify(f"{result.element}: {len(result.lines)} lines measured — synthesising…")
        self.run_synth(request, observation, (result.snr, result.cont_pct))

    # --- step 1: per-line synthesis (streamed) -----------------------------------------
    def run_synth(self, request=None, observation=None, uncertainty=None):
        res = self.selected()
        if res is None:
            self.notify("Select an abundance object first.")
            return
        self.setup.element = res.element
        ctx = request if request is not None else self.setup.context()
        if ctx is None:
            self.notify(self.setup.reason or "Configure the synthesis setup first.")
            return
        self._settings_uid = None
        self.push_setup_note()
        self.notify(f"Synthesising {res.element} lines…")
        obs = observation if observation is not None else self.get_observation(res.ref_uid)
        if obs is None:
            self.notify("The reference observation is no longer loaded.")
            return
        if request is None:
            self.draft.commit()
            pulse(self._modal, "close")
        if uncertainty is None:
            uncertainty = (res.snr, res.cont_pct)
        snr, cont_pct = uncertainty
        kwargs = {"target": obs, "measured": res, "method": "profile", "profile": res.profile,
                  "snr": snr, "cont_pct": cont_pct,
                  "background_lines": self.setup.background_lines}
        self.jobs.submit_operation(res.name, self.task_kind, ctx,
                                   self._on_synth_done, kwargs=kwargs, result=res)

    task_kind = "abundance.fit"
    import_task_result = lambda self, result, handle: self._on_synth_done(
        self._adopt(result, handle.label))

    def _on_synth_done(self, result):
        self.push_results()
        if self.selected() is result:
            self._render_selected()
        self.notify(f"{result.element}: abundance and uncertainties computed.")

    # --- step 2: uncertainties ---------------------------------------------------------
    def run_uncert(self, res=None):
        """The error budget for ``res`` (the selected result by default).

        The synthesis names its own result when it chains into this, so a run that finishes
        while another object is selected still gets its uncertainties. The SNR and continuum
        uncertainty are the ones stored on ``res`` when its windows were measured; the modal's
        fields only apply to the next measurement.
        """
        res = res if res is not None else self.selected()
        if res is None:
            self.notify("Select an abundance object first.")
            return
        if res.stage == "lines":
            self.notify("Synthesise the lines before computing uncertainties.")
            return
        spec = self.get_observation(getattr(res, "ref_uid", None))
        if spec is None:
            self.notify("The reference observation for this result is no longer loaded.")
            return
        ab.finalize_result(res, spec, res.snr, res.cont_pct)
        s = res.summary
        self.push_results()
        if self.selected() is res:
            self._render_selected()
        if np.isfinite(s["a_mean"]):
            self.notify(f"{res.element}: ⟨{self.value_label(res)}⟩ = "
                        f"{s['a_mean']:.3f} ± {s['sigma_mean']:.3f} "
                        f"(phot {s['sigma_phot']:.3f}, cont {s['sigma_cont']:.3f}, n={s['n']})")

    # --- rendering ----------------------------------------------------------------------
    def _render_scatter(self):
        """Repaint both diagnostic panels from the lines as they stand.

        Called on every streamed row too, so the panels fill in during the synthesis: a line
        enters as soon as it has an abundance, with a zero-length bar until the uncertainties
        are computed and its σ(A) is known.
        """
        res = self.selected()
        inc = [(i, ln) for i, ln in enumerate(res.lines)
               if ln.include and np.isfinite(ln.a)] if res else []
        if not inc:
            self.src.data = self.empty_points()
            self._set_mean(None)
            return
        for fig in self.figs:
            fig.yaxis.axis_label = self.value_label(res)
        data = self._pack_points(inc)
        self.src.data = data
        self._crop_panels((data,))
        s = res.summary
        self._set_mean(s.get("a_mean", np.nan), s.get("sigma_mean", 0.0))

    def ref_html(self):
        if self._settings_uid:
            res = self._find(self._settings_uid)
            if res is not None:
                return f'<div class="hint">Reusing measured widths from {html.escape(res.ref_name)}. Start from the line mask set to measure again.</div>'
        obs = self.list_observations()
        cur = self.ref_uid or (obs[0][0] if obs else None)
        opts = "".join(
            f'<option value="{html.escape(str(uid))}"'
            f'{" selected" if uid == cur else ""}>{html.escape(name)}</option>'
            for uid, name in obs) or '<option value="">— no observation loaded —</option>'
        return ('<div class="field-group"><label>Target spectrum</label>'
                f'<select class="field-input abund-ref-select" aria-label="Target spectrum">{opts}</select></div>'
                + self.ctx.canvas.preparation_html(cur))

    def results_html(self):
        if not self.results:
            return ""
        rows = self._results_rows("A")
        return ('<div class="section-header result-section-head">Abundances</div>'
                + "".join(rows))

    def table_html(self, res):
        head = self._table_head((
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), (LOCAL_CONTINUUM_HEAD, None), ("blend", "blend"),
            (html.escape(self.value_label(res)), "a"),
            ("vshift", "vshift"), ("vdop", "vdop"), ("flags", None)))
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
            if self.show_toggled_only and not ln.include:
                continue          # header mask: keep only toggled-in lines (index i unchanged)
            chk = " checked" if ln.include else ""
            cls = " row-selected" if self._jumped == i else ""
            cls += "" if ln.include else " line-excluded"
            ew = self._fmt(ln.ew, ln.ew_err, digits=1)   # EW error known from step 1
            blend = self._fmt(ln.blend_frac, digits=2)
            a = self._fmt(ln.a, ln.sigma_a if res.stage == "uncert" else None, plus=True)
            vshift = self._fmt(ln.vshift, digits=2)
            vdop = self._fmt(ln.vdop, digits=2)
            fl = ", ".join(ln.flags)
            body.append(
                f'<tr class="abund-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="line-include-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'{self.continuum_cell(ln, i)}<td>{blend}</td><td>{a}</td><td>{vshift}</td><td>{vdop}</td>'
                f'<td class="abund-flags">{html.escape(fl)}</td></tr>')
        foot = ""
        if res.stage == "uncert" and np.isfinite(res.summary.get("a_mean", np.nan)):
            s = res.summary
            stale = bool(self.stale_reasons(res))
            foot = (f'<div class="abund-summary{" stale-block" if stale else ""}">'
                    + ('<span class="stale-tag">Previous solution</span> ' if stale else "")
                    + f'⟨{html.escape(self.value_label(res))}⟩ = '
                    f'{s["a_mean"]:+.3f} ± {s["sigma_mean"]:.3f}'
                    f'<span class="abund-budget"> (photon {s["sigma_phot"]:.3f} · '
                    f'continuum {s["sigma_cont"]:.3f} · scatter {s["sigma_scatter"]:.3f} · '
                    f'n={s["n"]})</span></div>')
        return (f'<table class="linelist-table abund-table abund-lines" data-tool="abund">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body)}</tbody></table>{foot}')
