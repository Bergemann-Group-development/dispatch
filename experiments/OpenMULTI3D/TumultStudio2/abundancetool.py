"""The Abundance tool: a mask-free, uncertainty-aware **3-step wizard** for element abundances.

Toolbar **Abund** opens a right-sidebar panel whose controls follow a step switcher shown in a
footer inside the centre panel:

1. **Line selection** (stellar-parameter-free): pick a reference observation, an element and an
   EW profile, and measure every candidate line's EW + FWHM. The result appears as an object in
   the left sidebar and a split view (spectrum on top, per-line table below) in the centre. The
   user curates the lines — clicking a row jumps the figure to that line over a window set by a
   multiple of its fitted FWHM.
2. **Synthesis**: using the tool's own synthesis setup (a :class:`synthsetup.SynthSetup`), build
   a per-line curve of growth and invert each EW to an abundance — **line by line, streaming the
   results into the table** as each lands. Jumping a row now also overlays the line's best-fit
   synthetic on the observation.
3. **Uncertainties**: set SNR / continuum % / scatter x-axis and compute each line's EW and
   abundance error budget plus the dark-themed A(X) scatter (which replaces the spectrum pane).

The numerical engine lives in :mod:`abundance` (pure, unit-tested) and the machinery shared
with the corrections wizard in :mod:`abundwizard`; this module is the abundance-specific
orchestration that drives ``synth.compute`` on the setup's atmosphere.
"""
from __future__ import annotations

import html
import logging
import uuid

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, Whisker

import abundance as ab
import store
import synth
from abundwizard import (PROFILES, SYNTH_PAD_AA, X_AXES, AbundanceWizard,
                         measure_candidate_lines)
from canvas import ACCENT
from controls import Field, Form
from linenav import load_master

DEFAULTS = {
    "element": "Fe", "profile": "gauss", "ew_min": 5.0,
    "snr": 100.0, "cont_pct": 0.5, "xaxis": "evi",
}


class AbundanceTool(AbundanceWizard):
    name = "abundance"
    prefix = "abund"
    kind = "result"
    btn_selector = "#abundance-btn"
    empty_table_msg = "No abundance result selected"
    y_axis_label = "A(X)"
    final_stage = "uncert"

    def __init__(self, ctx, *, setup, list_observations, get_observation, master_path):
        self.list_observations = list_observations  # () -> [(uid, name)]
        self.get_observation = get_observation      # (uid) -> Spectrum (normalised) | None
        self.master_path = master_path

        saved = dict(store.get_section("abundance"))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # step 1: what to measure (independent of stellar parameters). The EW error inputs
        # (SNR / continuum-%) live here too, so each line's EW uncertainty is known at
        # selection; a line whose photon EW error exceeds its EW is flagged automatically.
        self.form1 = Form("abund1", [
            Field("element", "Element (species)", kind="text", value=vals["element"]),
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=PROFILES),
            Field("ew_min", "Min EW [mÅ]", kind="number", value=vals["ew_min"], min=0),
            Field("snr", "SNR (observation)", kind="number", value=vals["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=vals["cont_pct"], min=0, step=0.1),
        ], on_change=self._on_field)
        # step 3: the scatter axis (SNR / continuum-% now live in step 1)
        self.form3 = Form("abund3", [
            Field("xaxis", "Scatter x-axis", kind="select", value=vals["xaxis"],
                  options=list(X_AXES)),
        ], on_change=self._on_field)

        super().__init__(ctx, setup=setup)
        self._ref_bridge = ctx.bridges.html("abund-ref")   # step-1 reference selector
        self.ref_uid = saved.get("ref_uid")

        # scatter with per-line error bars (the corrections wizard has none)
        src = ColumnDataSource(dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[], inc=[]))
        fig = self._build_fig(src)
        fig.add_layout(Whisker(base="x", upper="hi", lower="lo", source=src,
                               line_color=ACCENT, line_width=1.2))

    # --- Tool protocol -----------------------------------------------------------------
    def _on_activate_extra(self):
        self._ref_bridge.value = self.ref_html()    # step-1 reference selector

    def refresh(self):
        """Keep the step-1 reference selector in sync as observations load/unload."""
        if self.active:
            self._ref_bridge.value = self.ref_html()

    def on_filter_changed(self):
        """Re-measure the selected result's step-1 EWs when the canvas Filtered switch toggles.

        ``get_observation`` already returns the filter-effective reference, so this just re-runs
        the (cheap) EW measurement over each line's stored fit window and repaints the table — the
        EW tool does the same live via :meth:`EquivalentWidth.update`. Restricted to the ``"lines"``
        stage on purpose: once synthesis has run, the abundances were inverted against the EW
        measured at selection, so silently re-measuring under a new filter would misrepresent them.
        Windows / flags / include-toggles are left untouched, so user curation survives."""
        if not self.active:
            return
        res = self.selected()
        if res is None or res.stage != "lines":
            return
        spec = self.get_observation(getattr(res, "ref_uid", None))
        if spec is None or spec.yunit != u.dimensionless_unscaled:
            return
        v = self.form1.values()
        snr = float(v["snr"]) or None
        cont_pct = float(v["cont_pct"])
        changed = False
        for ln in res.lines:
            if not np.isfinite(ln.fit_window[0]):
                continue          # a no-window line carries no EW to refresh
            try:
                prof = ab.measure_ew(spec, window=ln.fit_window, snr=snr,
                                     cont_pct=cont_pct, profile=res.profile, x0=ln.lam0)
            except Exception:
                continue
            ln.ew, ln.ew_err = prof.mAA, prof.width_err_mAA
            changed = True
        if changed:
            self._render_selected()

    # --- session persistence (whole computed results, saved into the session slot) ------
    def snapshot(self):
        """Serialise the computed results + selection for the session slot. Each result records
        its reference observation by name too, so it can be re-linked after restore mints new uids."""
        obs = dict(self.list_observations())    # uid -> name
        out = []
        for r in self.results:
            d = r.to_dict()
            d["ref_name"] = obs.get(getattr(r, "ref_uid", None))
            out.append(d)
        return {"results": out, "selected": self.selected_uid}

    def restore(self, data):
        """Rebuild the results (and left-sidebar list) from a saved snapshot."""
        data = data or {}
        by_name = {}
        for uid, name in self.list_observations():
            by_name.setdefault(name, uid)   # first match wins on duplicate names
        results = []
        for d in (data.get("results") or []):
            res = ab.AbundanceResult.from_dict(d)
            # observation uids are re-minted on restore, so re-link by the saved reference name
            new_uid = by_name.get(d.get("ref_name"))
            if new_uid is not None:
                res.ref_uid = new_uid
            results.append(res)
        self.results = results
        uids = {r.uid for r in results}
        self.selected_uid = data.get("selected") if data.get("selected") in uids else None
        self.push_results()

    # --- state persistence -------------------------------------------------------------
    def _on_field(self, key, value):
        store.update_pref("abundance", key, value)
        if key == "xaxis":
            self.fig.xaxis.axis_label = X_AXES[str(value)][0]
            self._render_scatter()  # re-place points on the new axis

    def _extra_action(self, kind, arg):
        if kind == "ref":
            self.ref_uid = arg or None
            store.update_pref("abundance", "ref_uid", self.ref_uid or "")
        elif kind == "run-lines":
            self.run_lines()
        elif kind == "run-synth":
            self.run_synth()
        elif kind == "run-uncert":
            self.run_uncert()

    # --- shared-machinery hooks -----------------------------------------------------------
    def _summarise(self, lines):
        return ab.aggregate(lines)

    def _result_label(self, res):
        s = res.summary
        return (f"{res.element}  ⟨A⟩={s['a_mean']:.2f}±{s['sigma_mean']:.2f}"
                if np.isfinite(s["a_mean"]) else f"{res.element} · uncert")

    def _jump_decor(self, res, ln):
        self.canvas.set_abund_overlay(ln.fit_spec, ln.window)
        self.canvas.set_abund_hatch(ln.fit_func, ln.window)
        # snap the window bounds between the reference's points so membership is unambiguous
        spec = self.get_observation(getattr(res, "ref_uid", None))
        xpts = np.asarray(spec.values[0], dtype=float) if spec is not None else None
        self.canvas.set_abund_windows(ln.fit_window, ln.window, xpts)

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

    @staticmethod
    def _species(raw):
        """Parse an element field: ("Fe II", "Fe", stage_named) — collapse inner whitespace."""
        element = " ".join(str(raw).split())
        symbol = element.split()[0] if element else ""
        return element, symbol, len(element.split()) > 1

    # --- step 1: diagnostic-line selection ---------------------------------------------
    def run_lines(self):
        ref_uid, ref_spec, err = self._ref()
        if err:
            self.notify(err)
            return
        v = self.form1.values()
        element, symbol, stage_named = self._species(v["element"])
        if not symbol:
            self.notify("Enter an element, e.g. 'Fe I'.")
            return
        profile = str(v["profile"])
        ew_min = float(v["ew_min"])
        snr = float(v["snr"]) or None
        cont_pct = float(v["cont_pct"])
        self._error = None
        self.notify(f"Measuring {element} lines on the reference…")
        loader = lambda token=None: self._measure_lines(
            ref_uid, ref_spec, element, symbol, stage_named, profile, ew_min, snr, cont_pct)
        self.ctx.submit(f"abund-lines:{element}", loader, self._on_lines_done)

    def _measure_lines(self, ref_uid, ref_spec, element, symbol, stage_named, profile, ew_min,
                       snr, cont_pct):
        try:
            lo, hi = self._ref_range(ref_spec)
            master = load_master(self.master_path)
            ion = master.ion.astype(str)
            match = (ion == element) if stage_named else (ion.str.split().str[0] == symbol)
            cand = master[match & (master.lam0 >= lo) & (master.lam0 <= hi)]
            cand = cand.sort_values("lam0")
            if cand.empty:
                self._error = f"No {element} lines in the reference range {lo:.1f}–{hi:.1f} Å."
                return None
            lines = measure_candidate_lines(ref_spec, cand, profile=profile, ew_min=ew_min,
                                            snr=snr, cont_pct=cont_pct)
            if not lines:
                self._error = f"No {element} lines above {ew_min:.0f} mÅ in the reference."
                return None
            res = ab.AbundanceResult(
                element=element, symbol=symbol, lines=lines, summary=ab.aggregate([]),
                stage="lines", profile=profile, uid=uuid.uuid4().hex[:8])
            res.ref_uid = ref_uid
            res.name = f"{element} · {len(lines)} lines"
            return res
        except Exception as exc:
            logging.exception("Abundance line selection failed")
            self._error = str(exc)
            return None

    def _on_lines_done(self, result):
        if result is None:
            if self._error:
                self.notify(f"Line selection failed: {self._error}")
            return
        self.results.append(result)
        self.selected_uid = result.uid
        self.set_step(1)
        self.select_result(result.uid)   # focus it: enters abund-mode, renders the split view
        self.notify(f"{result.element}: {len(result.lines)} candidate lines — review and curate.")

    # --- step 2: per-line synthesis (streamed) -----------------------------------------
    def run_synth(self):
        res = self.selected()
        if res is None:
            self.notify("Select an abundance object first.")
            return
        ctx = self.setup.context()
        if ctx is None:
            self.notify(self.setup.reason or "Configure the synthesis setup (step 2).")
            return
        self._error = None
        self.notify(f"Synthesising {res.element} lines…")
        loader = lambda token=None: self._synth_lines(res, ctx, token)
        self.ctx.submit(f"abund-synth:{res.element}", loader, self._on_synth_done)

    def _synth_lines(self, res, ctx, token):
        try:
            prefs, grid, base_rv = ctx["prefs"], ctx["grid"], dict(ctx["run_vals"])
            target, atmos_file = ctx["target"], ctx["atmos_file"]
            feh = target[2]
            relative = bool(base_rv.get("relative_abundances", False))
            solar = synth.abund_table(prefs)
            base_abund = str(base_rv.get("abundances", "")).strip()
            symbol = res.symbol
            obs = self.get_observation(getattr(res, "ref_uid", None))
            if obs is None:
                self._error = "The reference observation for this result is no longer loaded."
                return None

            def synth_window(rv_over, window):
                rv = dict(base_rv)
                rv.update(rv_over)
                return synth.compute(target, prefs, rv, window[0] - SYNTH_PAD_AA,
                                     window[1] + SYNTH_PAD_AA, grid=grid, atmos_file=atmos_file,
                                     cancel_token=token)

            def abund_terms(a):
                val = a if relative else synth.xfe_to_abund(a, symbol, feh, solar)
                return " ".join(t for t in (base_abund, f"{symbol}={val:.4f}") if t)

            for idx, ln in enumerate(res.lines):
                if token is not None and token.cancelled:
                    raise synth.CancelledError()
                if not ln.include:
                    continue          # skip step-1-rejected / user-excluded lines (keeps flags)

                # the query species in isolation (weak-line screen), the full synthesis (blend +
                # gradient), the element-removed window (blend), and the adaptive abundance fit.
                def synth_only_fn(a, _w=ln.window):
                    return synth_window({"abundances": abund_terms(a),
                                         "relative_abundances": relative,
                                         "only_elements": symbol}, _w)

                def synth_all_fn(a, _w=ln.window):
                    return synth_window({"abundances": abund_terms(a),
                                         "relative_abundances": relative}, _w)

                # the dREW/dA flanks (a±grad_step full syntheses) fanned across one pool
                def synth_all_batch(a_list, _w=ln.window):
                    rvs = [{**base_rv, "abundances": abund_terms(a),
                            "relative_abundances": relative} for a in a_list]
                    return synth.compute_batch(target, prefs, rvs, _w[0] - SYNTH_PAD_AA,
                                               _w[1] + SYNTH_PAD_AA, grid=grid,
                                               atmos_file=atmos_file, cancel_token=token)

                def synth_excl_fn(_w=ln.window):
                    return synth_window({"exclude_elements": symbol}, _w)

                def abund_fit_fn(guess, _w=ln.window, _fw=ln.fit_window):
                    return synth.fit_line(obs, target, prefs, base_rv, symbol, guess,
                                          _w[0] - SYNTH_PAD_AA, _w[1] + SYNTH_PAD_AA, _fw,
                                          grid=grid, atmos_file=atmos_file, cancel_token=token)

                try:
                    ab.fit_line(ln, synth_only_fn=synth_only_fn, synth_all_fn=synth_all_fn,
                                synth_excl_fn=synth_excl_fn, abund_fit_fn=abund_fit_fn,
                                ew_meas=ln.ew, a_guess=0.0, profile=res.profile,
                                grad_batch_fn=synth_all_batch)
                except synth.CancelledError:
                    raise
                except Exception:
                    logging.exception("Abundance: fit for line %.3f failed", ln.lam0)
                    ln.flags = ("synth_failed",)
                res.relative = relative
                res.teff, res.logg, res.feh, res.vmic = target[0], target[1], feh, target[3]
                self.ctx.schedule(lambda r=res: self._stream_row(r))  # refresh the table live
            return res
        except synth.CancelledError:
            raise
        except Exception as exc:
            logging.exception("Abundance synthesis failed")
            self._error = str(exc)
            return None

    def _on_synth_done(self, result):
        if result is None:
            if self._error:
                self.notify(f"Synthesis failed: {self._error}")
            return
        result.stage = "synth"
        result.name = f"{result.element} · {len(result.lines)} lines · synth"
        self.push_results()
        if self.selected() is result:
            self.push_table(self.table_html(result))
        self.notify(f"{result.element}: synthesis done — set SNR/continuum in step 3.")

    # --- step 3: uncertainties ---------------------------------------------------------
    def run_uncert(self):
        res = self.selected()
        if res is None:
            self.notify("Select an abundance object first.")
            return
        if res.stage == "lines":
            self.notify("Run the synthesis (step 2) before computing uncertainties.")
            return
        spec = self.get_observation(getattr(res, "ref_uid", None))
        if spec is None:
            self.notify("The reference observation for this result is no longer loaded.")
            return
        v = self.form1.values()          # SNR / continuum-% are set in step 1
        snr = float(v["snr"]) or None
        cont_pct = float(v["cont_pct"])
        for ln in res.lines:
            if np.isfinite(ln.drewda):   # only lines that passed screening + were fit in step 2
                ab.finalize_line(ln, spec, snr, cont_pct, profile=res.profile)
        res.snr, res.cont_pct = snr, cont_pct
        res.summary = ab.aggregate(res.lines)
        res.stage = "uncert"
        s = res.summary
        res.name = self._result_label(res)
        self.push_results()
        self._render_selected()
        if np.isfinite(s["a_mean"]):
            self.notify(f"{res.element}: ⟨A⟩ = {s['a_mean']:.3f} ± {s['sigma_mean']:.3f} "
                        f"(phot {s['sigma_phot']:.3f}, cont {s['sigma_cont']:.3f}, n={s['n']})")

    @staticmethod
    def _ref_range(ref_spec):
        x, _ = ref_spec.values
        return float(np.min(x)), float(np.max(x))

    # --- rendering ----------------------------------------------------------------------
    def _render_scatter(self):
        res = self.selected()
        inc = [ln for ln in res.lines if ln.include and np.isfinite(ln.sigma_a)] if res else []
        if not inc:
            self.src.data = dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[], inc=[])
            self.mean_line.visible = self.mean_band.visible = False
            return
        xf = X_AXES[str(self.form3.value("xaxis"))][1]
        xs = [xf(ln) for ln in inc]
        los = [ln.a - ln.sigma_a for ln in inc]
        his = [ln.a + ln.sigma_a for ln in inc]
        self.src.data = dict(
            x=xs, a=[ln.a for ln in inc], lo=los, hi=his,
            lam=[ln.lam0 for ln in inc], ion=[ln.ion for ln in inc],
            inc=[ln.include for ln in inc])
        self._crop_scatter(xs, los + his)
        s = res.summary
        if np.isfinite(s.get("a_mean", np.nan)):
            self.mean_line.location = s["a_mean"]
            self.mean_line.visible = True
            self.mean_band.top = s["a_mean"] + s["sigma_mean"]
            self.mean_band.bottom = s["a_mean"] - s["sigma_mean"]
            self.mean_band.visible = True
        else:
            self.mean_line.visible = self.mean_band.visible = False

    def ref_html(self):
        obs = self.list_observations()
        cur = self.ref_uid or (obs[0][0] if obs else None)
        opts = "".join(
            f'<option value="{html.escape(str(uid))}"'
            f'{" selected" if uid == cur else ""}>{html.escape(name)}</option>'
            for uid, name in obs) or '<option value="">— no observation loaded —</option>'
        return ('<div class="field-group"><label>Reference (observation)</label>'
                f'<select class="field-input abund-ref-select">{opts}</select></div>')

    def results_html(self):
        if not self.results:
            return ""
        rows = self._results_rows("A", "", "data-result", "data-result-remove")
        return ('<div class="section-header result-section-head">Abundances</div>'
                + "".join(rows))

    def table_html(self, res):
        mask = " checked" if self.show_toggled_only else ""
        head = ('<tr><th><input type="checkbox" class="abund-mask-chk"'
                f'{mask} title="Show only included lines"></th>'
                "<th>λ₀ [Å]</th><th>ion</th><th>E<sub>low</sub></th><th>log gf</th>"
                "<th>EW [mÅ]</th><th>blend</th><th>A(X)</th><th>vshift</th><th>vdop</th>"
                "<th>flags</th></tr>")
        body = []
        for i, ln in enumerate(res.lines):
            if self.show_toggled_only and not ln.include:
                continue          # header mask: keep only toggled-in lines (index i unchanged)
            chk = " checked" if ln.include else ""
            cls = "" if ln.include else " line-excluded"
            ew = self._fmt(ln.ew, ln.ew_err, digits=1)   # EW error known from step 1
            blend = self._fmt(ln.blend_frac, digits=2)
            a = self._fmt(ln.a, ln.sigma_a if res.stage == "uncert" else None, plus=True)
            vshift = self._fmt(ln.vshift, digits=2)
            vdop = self._fmt(ln.vdop, digits=2)
            fl = ", ".join(ln.flags)
            body.append(
                f'<tr class="abund-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="abund-line-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'<td>{blend}</td><td>{a}</td><td>{vshift}</td><td>{vdop}</td>'
                f'<td class="abund-flags">{html.escape(fl)}</td></tr>')
        foot = ""
        if res.stage == "uncert" and np.isfinite(res.summary.get("a_mean", np.nan)):
            s = res.summary
            foot = (f'<div class="abund-summary">⟨A({html.escape(res.symbol)})⟩ = '
                    f'{s["a_mean"]:+.3f} ± {s["sigma_mean"]:.3f}'
                    f'<span class="abund-budget"> (photon {s["sigma_phot"]:.3f} · '
                    f'continuum {s["sigma_cont"]:.3f} · scatter {s["sigma_scatter"]:.3f} · '
                    f'n={s["n"]})</span></div>')
        return (f'<table class="linelist-table abund-table"><thead>{head}</thead>'
                f'<tbody>{"".join(body)}</tbody></table>{foot}')
