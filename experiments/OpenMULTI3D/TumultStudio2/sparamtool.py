"""The Stellar-Parameter tool: a 3-step wizard determining (Teff, logg, [Fe/H], vmic)
from a target observation via Fe I/II excitation–ionization balance.

Toolbar **Params** opens a right-sidebar panel whose controls follow a step switcher shown
in a footer inside the centre panel (the third sibling of the Abund / Corr wizards, on the
same :class:`abundwizard.AbundanceWizard` machinery):

1. **Fe line selection** (stellar-parameter-free): pick a reference observation and an EW
   profile, and measure every Fe I + Fe II candidate line's EW + FWHM
   (:func:`abundwizard.measure_candidate_lines` — the Abund Tool's measuring loop with the
   element fixed to both Fe stages). The user curates the lines in the split view.
2. **Solve.** Using the tool's own :class:`synthsetup.SynthSetup` (whose Teff/logg/[Fe/H]/
   vmic form is the *initial guess*), iterate: at the current parameters invert each line's
   EW to A(Fe) through a species-only curve of growth (``synth.compute_batch`` ladder →
   :func:`abundance.build_cog` → :func:`abundance.cog_invert`, exactly the Corr tool's
   inner loop), compute the four balance diagnostics, and take a damped quasi-Newton step
   (:mod:`sparams`) until all four are inside tolerance. Rows and an iteration status line
   stream into the table as they land. The first iteration also blend-screens each line
   (full vs Fe-removed synthesis, as the Abund tool does).
3. **Diagnostics.** The per-line A(Fe) scatter (Fe I accent, Fe II red) over E_low /
   reduced EW / λ with the weighted regression line, plus the converged parameters ± σ.
"""
from __future__ import annotations

import html
import logging
import uuid

import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, Slope, Whisker

import abundance as ab
import sparams as sp
import store
import synth
from abundwizard import (PROFILES, SYNTH_PAD_AA, X_AXES, AbundanceWizard,
                         measure_candidate_lines)
from canvas import ACCENT
from controls import Field, Form
from linenav import load_master
from tumult.abundance import xfe_to_abund

FE2_COLOR = "#e05c5c"      # Fe II scatter points (Fe I uses the shared accent)
BLEND_MAX = 0.30           # iteration-1 blend screen: contaminant EW fraction ceiling
VMIC_BOUNDS = (0.1, 5.0)   # vmic is a synthesis input, not a grid axis — keep it sane


def _is_shutdown(exc):
    """True for a ProcessPoolExecutor refusing new work because the interpreter is tearing
    down (``cannot schedule new futures after interpreter shutdown``) — the background solve
    thread outliving its Bokeh session. Treated as a cancel, not a per-line synthesis failure,
    so the loop aborts instead of re-issuing (and re-failing) every remaining line."""
    return isinstance(exc, RuntimeError) and "interpreter shutdown" in str(exc)

DEFAULTS = {
    "profile": "gauss", "ew_min": 5.0, "snr": 100.0, "cont_pct": 0.5, "xaxis": "evi",
    "max_iter": 8.0, "span": 0.4, "points": 3.0,
    "tol_exc": sp.DEFAULT_TOLS["s_exc"], "tol_rew": sp.DEFAULT_TOLS["s_rew"],
    "tol_ion": sp.DEFAULT_TOLS["d_ion"], "tol_feh": sp.DEFAULT_TOLS["d_feh"],
}


class SParamTool(AbundanceWizard):
    name = "sparams"
    prefix = "sparam"
    kind = "sparam"
    btn_selector = "#sparams-btn"
    section = "sparams"
    empty_table_msg = "No parameter result selected"
    y_axis_label = "A(Fe)"
    final_stage = "solved"
    extra_sort_keys = {"blend": lambda ln: ln.blend_frac}
    extra_sort_asc = {"blend": True}

    def __init__(self, ctx, *, setup, list_observations, get_observation, master_path):
        self.list_observations = list_observations  # () -> [(uid, name)]
        self.get_observation = get_observation      # (uid) -> Spectrum (normalised) | None
        self.master_path = master_path

        saved = dict(store.get_section(self.section))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # step 1: the measurement (element fixed: Fe I + Fe II)
        self.form1 = Form("sparam1", [
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=PROFILES),
            Field("ew_min", "Min EW [mÅ]", kind="number", value=vals["ew_min"], min=0),
            Field("snr", "SNR (observation)", kind="number", value=vals["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=vals["cont_pct"], min=0, step=0.1),
        ], on_change=self._on_field)
        # step 2: the solver knobs (the initial guess lives in the SynthSetup run form)
        self.form2 = Form("sparam2", [
            Field("max_iter", "Max iterations", kind="number", value=vals["max_iter"],
                  min=1, step=1),
            Field("span", "CoG half-span [dex]", kind="number", value=vals["span"],
                  min=0.05, step=0.05),
            Field("points", "CoG points", kind="number", value=vals["points"],
                  min=2, step=1),
            Field("tol_exc", "Tol: excitation slope [dex/eV]", kind="number",
                  value=vals["tol_exc"], min=0, step=0.001),
            Field("tol_rew", "Tol: REW slope [dex/dex]", kind="number",
                  value=vals["tol_rew"], min=0, step=0.005),
            Field("tol_ion", "Tol: ionization ΔA [dex]", kind="number",
                  value=vals["tol_ion"], min=0, step=0.005),
            Field("tol_feh", "Tol: [Fe/H] consistency [dex]", kind="number",
                  value=vals["tol_feh"], min=0, step=0.005),
        ], on_change=self._on_field)
        # step 3: the scatter axis
        self.form3 = Form("sparam3", [
            Field("xaxis", "Scatter x-axis", kind="select", value=vals["xaxis"],
                  options=list(X_AXES)),
        ], on_change=self._on_field)

        super().__init__(ctx, setup=setup)
        self._ref_bridge = ctx.bridges.html("sparam-ref")   # step-1 reference selector
        self.ref_uid = saved.get("ref_uid")

        # Fe I points on the base figure's source; Fe II on a second source in red;
        # whiskers on both; a fitted regression line over the *Fe I* points.
        src = ColumnDataSource(dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[]))
        fig = self._build_fig(src)
        fig.add_layout(Whisker(base="x", upper="hi", lower="lo", source=src,
                               line_color=ACCENT, line_width=1.2))
        self.src2 = ColumnDataSource(dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[]))
        fig.scatter("x", "a", source=self.src2, size=8, fill_color=FE2_COLOR,
                    line_color="#fff", line_width=1)
        fig.add_layout(Whisker(base="x", upper="hi", lower="lo", source=self.src2,
                               line_color=FE2_COLOR, line_width=1.2))
        self.trend = Slope(gradient=0, y_intercept=0, line_color=ACCENT,
                           line_dash="dotted", line_width=1.4)
        self.trend.visible = False
        fig.add_layout(self.trend)

    # --- Tool protocol -----------------------------------------------------------------
    @property
    def models(self):
        return [*self.form1.models, *self.form2.models, *self.form3.models,
                self.actions, *self.setup.models]

    def _on_activate_extra(self):
        self._ref_bridge.value = self.ref_html()
        self.form2.render()

    def refresh(self):
        """Keep the step-1 reference selector in sync as observations load/unload."""
        if self.active:
            self._ref_bridge.value = self.ref_html()

    def on_filter_changed(self):
        """Re-measure the selected result's step-1 EWs when the canvas Filtered switch
        toggles — same rules as the Abund tool: only while still at the ``"lines"`` stage
        (a solved result's abundances were inverted against the EWs measured then)."""
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
                continue
            try:
                prof = ab.measure_ew(spec, window=ln.fit_window, snr=snr,
                                     cont_pct=cont_pct, profile=res.profile, x0=ln.lam0)
            except Exception:
                continue
            ln.ew, ln.ew_err = prof.mAA, prof.width_err_mAA
            changed = True
        if changed:
            self._render_selected()

    # --- session persistence -------------------------------------------------------------
    def snapshot(self):
        obs = dict(self.list_observations())
        out = []
        for r in self.results:
            d = r.to_dict()
            d["ref_name"] = obs.get(getattr(r, "ref_uid", None))
            out.append(d)
        return {"results": out, "selected": self.selected_uid}

    def restore(self, data):
        data = data or {}
        by_name = {}
        for uid, name in self.list_observations():
            by_name.setdefault(name, uid)
        results = []
        for d in (data.get("results") or []):
            res = sp.SParamResult.from_dict(d)
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
        store.update_pref(self.section, key, value)
        if key == "xaxis":
            self.fig.xaxis.axis_label = X_AXES[str(value)][0]
            self._render_scatter()

    def _extra_action(self, kind, arg):
        if kind == "ref":
            self.ref_uid = arg or None
            store.update_pref(self.section, "ref_uid", self.ref_uid or "")
        elif kind == "run-lines":
            self.run_lines()
        elif kind == "run-solve":
            self.run_solve()

    # --- shared-machinery hooks -----------------------------------------------------------
    def _toggle_line(self, token):
        """Include-toggle a line; on a solved result, re-aggregate the diagnostics and
        uncertainties from the stored final-iteration abundances (no re-synthesis)."""
        res = self.selected()
        if res is None:
            return
        try:
            idx = int(token)
        except ValueError:
            return
        if not (0 <= idx < len(res.lines)):
            return
        res.lines[idx].include = not res.lines[idx].include
        if res.stage == self.final_stage:
            self._reaggregate(res)
            res.name = self._result_label(res)
            self.push_results()
        self._render_selected()

    def _reaggregate(self, res):
        feh = float(res.params.get("feh", float("nan")))
        res.diags = sp.diagnostics(res.lines, feh, res.solar_fe)
        res.sigmas = sp.uncertainties(res.diags, res.jac or dict(sp.DEFAULT_JAC))

    def _result_label(self, res):
        p = res.params
        if res.stage == "solved" and p:
            mark = "" if res.converged else " (max iter)"
            return (f"Teff {p['teff']:.0f} · logg {p['logg']:.2f} · "
                    f"[Fe/H] {p['feh']:+.2f} · ξ {p['vmic']:.2f}{mark}")
        return res.name

    def _jump_decor(self, res, ln):
        # no per-line synthetic fit is kept during the solve, so decorate with the
        # step-1 profile hatch + windows (and the overlay only when one exists)
        if ln.fit_spec is not None:
            self.canvas.set_abund_overlay(ln.fit_spec, ln.window)
        if ln.fit_func is not None:
            self.canvas.set_abund_hatch(ln.fit_func, ln.window)
        spec = self.get_observation(getattr(res, "ref_uid", None))
        xpts = np.asarray(spec.values[0], dtype=float) if spec is not None else None
        self.canvas.set_abund_windows(ln.fit_window, ln.window, xpts)

    # --- reference resolution ----------------------------------------------------------
    def _ref(self):
        obs = self.list_observations()
        ref_uid = self.ref_uid or (obs[0][0] if obs else None)
        spec = self.get_observation(ref_uid) if ref_uid else None
        if spec is None:
            return None, None, "Load and select a target observation to measure EWs on."
        if spec.yunit != u.dimensionless_unscaled:
            return None, None, ("Normalise the target first — the Fe balance needs "
                                "continuum-normalised flux.")
        return ref_uid, spec, None

    # --- step 1: Fe I + Fe II line selection ---------------------------------------------
    def run_lines(self):
        ref_uid, ref_spec, err = self._ref()
        if err:
            self.notify(err)
            return
        v = self.form1.values()
        profile = str(v["profile"])
        ew_min = float(v["ew_min"])
        snr = float(v["snr"]) or None
        cont_pct = float(v["cont_pct"])
        self._error = None
        self.notify("Measuring Fe I + Fe II lines on the target…")
        loader = lambda token=None: self._measure_lines(
            ref_uid, ref_spec, profile, ew_min, snr, cont_pct)
        self.ctx.submit("sparam-lines:Fe", loader, self._on_lines_done)

    def _measure_lines(self, ref_uid, ref_spec, profile, ew_min, snr, cont_pct):
        try:
            x, _ = ref_spec.values
            lo, hi = float(np.min(x)), float(np.max(x))
            master = load_master(self.master_path)
            ion = master.ion.astype(str)
            match = ion.str.split().str[0] == "Fe"
            cand = master[match & (master.lam0 >= lo) & (master.lam0 <= hi)]
            cand = cand.sort_values("lam0")
            if cand.empty:
                self._error = f"No Fe lines in the target range {lo:.1f}–{hi:.1f} Å."
                return None
            lines = measure_candidate_lines(ref_spec, cand, profile=profile, ew_min=ew_min,
                                            snr=snr, cont_pct=cont_pct)
            if not lines:
                self._error = f"No Fe lines above {ew_min:.0f} mÅ in the target."
                return None
            res = sp.SParamResult(lines=lines, stage="lines", profile=profile,
                                  uid=uuid.uuid4().hex[:8])
            res.ref_uid = ref_uid
            n1 = sum(1 for ln in lines if sp.is_ion_stage(ln.ion, 1))
            n2 = sum(1 for ln in lines if sp.is_ion_stage(ln.ion, 2))
            res.name = f"Fe balance · {n1} Fe I / {n2} Fe II"
            return res
        except Exception as exc:
            logging.exception("Stellar-parameter line selection failed")
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
        self.select_result(result.uid)
        self.notify(f"{len(result.lines)} candidate Fe lines — curate them, then solve "
                    "in step 2 (a healthy solve wants Fe II lines too).")

    # --- step 2: the balance solve --------------------------------------------------------
    def run_solve(self):
        res = self.selected()
        if res is None:
            self.notify("Select a parameter result first (run step 1).")
            return
        ctx = self.setup.context()
        if ctx is None:
            self.notify(self.setup.reason or "Configure the synthesis setup (step 2).")
            return
        if self.get_observation(getattr(res, "ref_uid", None)) is None:
            self.notify("The target observation for this result is no longer loaded.")
            return
        self._error = None
        self.notify("Solving stellar parameters from the Fe balance…")
        loader = lambda token=None: self._solve(res, ctx, token)
        self.ctx.submit("sparam-solve:Fe", loader, self._on_solve_done)

    @staticmethod
    def _strip_fe(text):
        """Drop any user Fe= term from the base abundance overrides (the solver owns Fe)."""
        return " ".join(t for t in str(text or "").split()
                        if not t.lower().startswith("fe="))

    def _grid_bounds(self, grid, atmos_file):
        """Box bounds of the atmosphere grid (a hull proxy: the interpolation itself
        raises on a truly unreachable query, which the solve loop surfaces)."""
        if atmos_file or not grid.get("folder"):
            return {"vmic": VMIC_BOUNDS}
        try:
            files = synth.discover_atmos_grid(grid["folder"],
                                              grid.get("pattern") or synth.DEFAULT_PATTERN)
        except Exception:
            return {"vmic": VMIC_BOUNDS}
        keys = np.array([k[:3] for k in files], dtype=float)
        bounds = {name: (float(keys[:, i].min()), float(keys[:, i].max()))
                  for i, name in enumerate(("teff", "logg", "feh"))}
        vmics = [k[3] for k in files if len(k) > 3]
        bounds["vmic"] = ((float(min(vmics)), float(max(vmics)))
                          if vmics else VMIC_BOUNDS)
        return bounds

    def _solve(self, res, ctx, token):
        try:
            prefs, grid = ctx["prefs"], ctx["grid"]
            atmos_file = ctx["atmos_file"]
            base_rv = dict(ctx["run_vals"])
            relative = bool(base_rv.get("relative_abundances", False))
            base_abund = self._strip_fe(base_rv.get("abundances", ""))
            solar = synth.abund_table(prefs)
            solar_fe = solar.get("Fe")
            if solar_fe is None:
                self._error = ("No solar Fe abundance — configure the abundance file in "
                               "Preferences → Synthesis.")
                return None

            v2 = self.form2.values()
            max_iter = max(1, int(float(v2["max_iter"])))
            span = float(v2["span"])
            npts = max(2, int(float(v2["points"])))
            tols = {"s_exc": float(v2["tol_exc"]), "s_rew": float(v2["tol_rew"]),
                    "d_ion": float(v2["tol_ion"]), "d_feh": float(v2["tol_feh"])}

            params = {k: float(base_rv[k]) for k in sp.PARAM_KEYS}
            state = sp.SolverState(params=dict(params),
                                   bounds=self._grid_bounds(grid, atmos_file))
            res.params0 = dict(params)
            res.solar_fe = float(solar_fe)
            res.history = []
            res.converged = False

            def rv_point(a_rel, feh, only=True, exclude=False):
                rv = dict(base_rv)
                rv.pop("only_elements", None)
                rv.pop("exclude_elements", None)
                if exclude:
                    rv["abundances"] = base_abund
                    rv["exclude_elements"] = "Fe"
                else:
                    val = a_rel if relative else xfe_to_abund(a_rel, "Fe", feh, solar)
                    rv["abundances"] = " ".join(
                        t for t in (base_abund, f"Fe={val:.4f}") if t)
                    if only:
                        rv["only_elements"] = "Fe"
                rv["relative_abundances"] = relative
                return rv

            # each line's ladder centre carries over between iterations ([Fe/H]-relative)
            centres = {}
            diags = {}
            for it in range(1, max_iter + 1):
                if token is not None and token.cancelled:
                    raise synth.CancelledError()
                target = tuple(params[k] for k in sp.PARAM_KEYS)
                feh = params["feh"]
                for i, ln in enumerate(res.lines):
                    if token is not None and token.cancelled:
                        raise synth.CancelledError()
                    if not ln.include or not np.isfinite(ln.ew) or ln.ew <= 0:
                        continue
                    c = centres.get(i, 0.0)
                    ladder = np.linspace(c - span, c + span, npts)
                    screen = it == 1        # blend screen rides along on iteration 1
                    rv_list = [rv_point(a, feh) for a in ladder]
                    if screen:
                        rv_list += [rv_point(c, feh, only=False),   # full synthesis
                                    rv_point(c, feh, exclude=True)]  # Fe removed
                    try:
                        specs = synth.compute_batch(
                            target, prefs, rv_list,
                            ln.window[0] - SYNTH_PAD_AA, ln.window[1] + SYNTH_PAD_AA,
                            grid=grid, atmos_file=atmos_file, cancel_token=token)
                        if screen:
                            blend, _, _ = ab.blend_fraction(
                                specs[-2], specs[-1], window=ln.fit_window,
                                profile=res.profile, x0=ln.lam0)
                            ln.blend_frac = float(blend)
                            if np.isfinite(blend) and blend > BLEND_MAX:
                                if "blended" not in ln.flags:
                                    ln.flags = ln.flags + ("blended",)
                                ln.include = False
                                self.ctx.schedule(lambda r=res: self._stream_row(r))
                                continue
                            specs = specs[:-2]
                        a_grid, ew_grid = ab.build_cog(
                            None, ladder, window=ln.fit_window, profile=res.profile,
                            x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs: _s)
                        ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                        inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
                        if inv["out_of_range"]:
                            # recentre a doubled ladder on the extrapolated abundance and
                            # try once more — the guess may simply be far off at first
                            c2 = float(np.clip(inv["a"], c - 3 * span, c + 3 * span))
                            ladder2 = np.linspace(c2 - 2 * span, c2 + 2 * span, npts)
                            specs2 = synth.compute_batch(
                                target, prefs, [rv_point(a, feh) for a in ladder2],
                                ln.window[0] - SYNTH_PAD_AA, ln.window[1] + SYNTH_PAD_AA,
                                grid=grid, atmos_file=atmos_file, cancel_token=token)
                            a_grid, ew_grid = ab.build_cog(
                                None, ladder2, window=ln.fit_window, profile=res.profile,
                                x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs2: _s)
                            inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
                        if inv["out_of_range"]:
                            ln.a = float("nan")
                            ln.flags = tuple(f for f in ln.flags if f != "out_of_range") \
                                + ("out_of_range",)
                        else:
                            a_rel = float(inv["a"])
                            centres[i] = a_rel
                            ln.a = xfe_to_abund(a_rel, "Fe", feh, solar)
                            ln.sigma_a = float(inv["sigma_a"])
                            ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
                            ln.a_grid, ln.ew_grid = tuple(a_grid), tuple(ew_grid)
                            ln.method = "cog"
                            ln.flags = tuple(f for f in ln.flags if f != "out_of_range")
                    except synth.CancelledError:
                        raise
                    except Exception as exc:
                        # a cancel can surface as a non-CancelledError — a killed pool
                        # (BrokenProcessPool/RuntimeError) mid-submit, or the executor
                        # refusing work once the interpreter is exiting. Abort the whole
                        # solve then, rather than logging (and re-issuing) every remaining line.
                        if (token is not None and token.cancelled) or _is_shutdown(exc):
                            raise synth.CancelledError()
                        logging.exception("Fe balance: line %.3f failed", ln.lam0)
                        ln.a = float("nan")
                        ln.flags = ("synth_failed",)
                    self.ctx.schedule(lambda r=res: self._stream_row(r))

                diags = sp.diagnostics(res.lines, feh, solar_fe)
                res.history.append({"params": dict(params), "diags": dict(diags)})
                res.params = dict(params)
                res.diags = dict(diags)
                self.ctx.schedule(lambda r=res: self._stream_row(r))  # show the status line
                if sp.converged(diags, tols):
                    res.converged = True
                    break
                params = sp.step(state, diags)

            res.jac = dict(state.jac)
            res.sigmas = sp.uncertainties(diags, state.jac)
            res.snr = float(self.form1.value("snr")) or None
            res.cont_pct = float(self.form1.value("cont_pct"))
            return res
        except synth.CancelledError:
            raise
        except Exception as exc:
            logging.exception("Stellar-parameter solve failed")
            self._error = str(exc)
            return None

    def _on_solve_done(self, result):
        if result is None:
            if self._error:
                self.notify(f"Solve failed: {self._error}")
            return
        result.stage = "solved"
        result.name = self._result_label(result)
        self.push_results()
        if self.selected() is result:
            self._render_selected()
        p, s = result.params, result.sigmas
        state = "converged" if result.converged else "hit max iterations"
        self.notify(f"Fe balance {state} after {len(result.history)} it: "
                    f"Teff={p['teff']:.0f}±{s.get('teff', float('nan')):.0f} "
                    f"logg={p['logg']:.2f}±{s.get('logg', float('nan')):.2f} "
                    f"[Fe/H]={p['feh']:+.2f}±{s.get('feh', float('nan')):.2f} "
                    f"vmic={p['vmic']:.2f}±{s.get('vmic', float('nan')):.2f}",
                    level="success" if result.converged else "warn")

    # --- rendering ----------------------------------------------------------------------
    def _render_scatter(self):
        res = self.selected()
        inc = [ln for ln in res.lines if ln.include and np.isfinite(ln.a)] if res else []
        empty = dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[])
        if not inc:
            self.src.data = dict(empty)
            self.src2.data = dict(empty)
            self.mean_line.visible = self.mean_band.visible = False
            self.trend.visible = False
            return
        xf = X_AXES[str(self.form3.value("xaxis"))][1]

        def pack(lns):
            xs = [xf(ln) for ln in lns]
            sig = [ln.sigma_a if np.isfinite(ln.sigma_a) else 0.0 for ln in lns]
            return dict(x=xs, a=[ln.a for ln in lns],
                        lo=[ln.a - s for ln, s in zip(lns, sig)],
                        hi=[ln.a + s for ln, s in zip(lns, sig)],
                        lam=[ln.lam0 for ln in lns], ion=[ln.ion for ln in lns])

        fe1 = [ln for ln in inc if sp.is_ion_stage(ln.ion, 1)]
        fe2 = [ln for ln in inc if sp.is_ion_stage(ln.ion, 2)]
        d1, d2 = pack(fe1), pack(fe2)
        self.src.data = d1
        self.src2.data = d2
        self._crop_scatter(d1["x"] + d2["x"], d1["lo"] + d1["hi"] + d2["lo"] + d2["hi"])

        # weighted regression over the Fe I points on the chosen axis
        slope, _, intercept = sp.weighted_slope(
            d1["x"], d1["a"], [ln.sigma_a for ln in fe1])
        if np.isfinite(slope):
            self.trend.gradient, self.trend.y_intercept = slope, intercept
            self.trend.visible = True
        else:
            self.trend.visible = False

        diags = res.diags or {}
        a_mean = diags.get("a_mean", float("nan"))
        if np.isfinite(a_mean):
            err = diags.get("d_feh_err", float("nan"))
            err = err if np.isfinite(err) else 0.0
            self.mean_line.location = a_mean
            self.mean_line.visible = True
            self.mean_band.top = a_mean + err
            self.mean_band.bottom = a_mean - err
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
        return ('<div class="field-group"><label>Target (observation)</label>'
                f'<select class="field-input sparam-ref-select">{opts}</select></div>')

    def results_html(self):
        if not self.results:
            return ""
        rows = self._results_rows("P", " sparam-result-item", "data-sparam-result",
                                  "data-sparam-remove")
        return ('<div class="section-header result-section-head">Stellar parameters</div>'
                + "".join(rows))

    @staticmethod
    def _diag_cell(diags, key, tol, digits=3):
        v = diags.get(key, float("nan"))
        if not np.isfinite(v):
            return "—"
        mark = "✓" if abs(v) <= tol else "✗"
        return f"{v:+.{digits}f} {mark}"

    def _status_html(self, res):
        """The iteration/convergence block under the table (streams during the solve)."""
        if not res.history:
            return ""
        last = res.history[-1]
        p, d = last["params"], last["diags"]
        v2 = self.form2.values()
        tols = {"s_exc": float(v2["tol_exc"]), "s_rew": float(v2["tol_rew"]),
                "d_ion": float(v2["tol_ion"]), "d_feh": float(v2["tol_feh"])}
        head = (f'iter {len(res.history)} · Teff {p["teff"]:.0f} · logg {p["logg"]:.2f} · '
                f'[Fe/H] {p["feh"]:+.2f} · ξ {p["vmic"]:.2f}')
        cells = (f'E_low slope {self._diag_cell(d, "s_exc", tols["s_exc"])} · '
                 f'REW slope {self._diag_cell(d, "s_rew", tols["s_rew"])} · '
                 f'Fe II−I {self._diag_cell(d, "d_ion", tols["d_ion"])} · '
                 f'[Fe/H] cons. {self._diag_cell(d, "d_feh", tols["d_feh"])}'
                 f' (n = {d.get("n_fe1", 0)} + {d.get("n_fe2", 0)})')
        if res.stage == "solved":
            s = res.sigmas or {}
            fmt = self._fmt
            state = "converged" if res.converged else "max iterations reached"
            head = (f'<b>{state}</b> after {len(res.history)} it · '
                    f'Teff {fmt(res.params.get("teff"), s.get("teff"), digits=0)} K · '
                    f'logg {fmt(res.params.get("logg"), s.get("logg"), digits=2)} · '
                    f'[Fe/H] {fmt(res.params.get("feh"), s.get("feh"), plus=True, digits=2)} · '
                    f'ξ {fmt(res.params.get("vmic"), s.get("vmic"), digits=2)} km/s')
        return (f'<div class="abund-summary sparam-status">{head}'
                f'<span class="abund-budget"> — {cells}</span></div>')

    def table_html(self, res):
        head = self._table_head("sparam-mask-chk", (
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), ("blend", "blend"),
            ("A(Fe)", "a"), ("flags", None)))
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
            if self.show_toggled_only and not ln.include:
                continue
            chk = " checked" if ln.include else ""
            cls = "" if ln.include else " line-excluded"
            ew = self._fmt(ln.ew, ln.ew_err, digits=1)
            blend = self._fmt(ln.blend_frac, digits=2)
            a = self._fmt(ln.a, ln.sigma_a if res.stage == "solved" else None, plus=True)
            fl = ", ".join(ln.flags)
            body.append(
                f'<tr class="sparam-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="sparam-line-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'<td>{blend}</td><td>{a}</td>'
                f'<td class="abund-flags">{html.escape(fl)}</td></tr>')
        return (f'<table class="linelist-table abund-table sparam-lines">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body)}</tbody></table>{self._status_html(res)}')
