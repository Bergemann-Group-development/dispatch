"""The Abundance-Correction tool: a 3-step wizard for **NLTE abundance corrections vs 1D LTE**.

Toolbar **Corr** opens a right-sidebar panel whose controls follow a step switcher shown in a
footer inside the centre panel. It reuses the abundance pipeline numerics
(:mod:`abundance` / :mod:`tumult.abundance`), the :class:`synthsetup.SynthSetup` config cluster,
the machinery shared with the abundance wizard (:mod:`abundwizard`), and the app-global grid
registry — one grid of which is flagged as the **1D LTE reference**.

1. **Reference & lines.** Pick stellar parameters (or a direct atmosphere file — the NLTE target)
   and a single model atom. "Compute reference" interpolates the flagged 1D LTE grid at those
   parameters and synthesises the reference spectrum **straight from the model atom** (no line
   list, via ``synth.compute_reference``). That run's ``linelist_out.txt`` supplies the candidate
   lines; each line's EW/FWHM/window is measured on the reference spectrum
   (:func:`abundance.measure_line`).
2. **NLTE corrections.** For every line, build an NLTE curve of growth on the *picked* atmosphere
   (model atom only, departures + spectrum, ``synth.compute``) and invert it
   (:func:`tumult.abundance.cog_invert`) to the abundance that reproduces the reference EW. The
   reference sits at ``[X/Fe]=0``, so that abundance **is** the correction ΔA = A(NLTE) − A(LTE).
   Rows stream into the table as each lands.
3. **Scatter.** Plot the line-to-line scatter of the corrections (mean ± population σ).
"""
from __future__ import annotations

import html
import logging
import uuid
from pathlib import Path

import numpy as np
from bokeh.models import ColumnDataSource

import abundance as ab
import store
import synth
from abundwizard import PROFILES, SYNTH_PAD_AA, X_AXES, AbundanceWizard
from controls import Field, Form

DEFAULTS = {"profile": "gauss", "ew_min": 5.0, "xaxis": "evi"}

CORR_SPAN = 0.6        # curve-of-growth spans ±CORR_SPAN dex in [X/Fe] about the reference
CORR_N = 7             # number of ladder points (kept small — each is an NLTE solve)
MIN_FIT_POINTS = 7     # a fit region sampled by fewer reference points than this is too sparse
                       # to trust — flag "few_points" and exclude the line by default

# Each correction carries a formal error after step 2 — the reference EW's
# fit-covariance error (model mismatch, nonzero even on a noiseless synthetic spectrum)
# propagated through the CoG slope — so the element mean is the inverse-variance
# weighted `aggregate`. `aggregate_unweighted` remains the fallback for results without
# per-line errors (step 1 only, or restored from a pre-error session).
from tumult.abundance import aggregate, aggregate_unweighted  # noqa: E402


def _summary(lines):
    weighted = aggregate(lines)
    if weighted["n"] > 0:
        return weighted
    return aggregate_unweighted(lines)


class AbundCorrTool(AbundanceWizard):
    name = "abundcorr"
    prefix = "corr"
    kind = "corr"
    btn_selector = "#abundcorr-btn"
    empty_table_msg = "No correction selected"
    y_axis_label = "ΔA (NLTE − LTE)"
    final_stage = "corr"

    def __init__(self, ctx, *, setup, grids, atom_provider, add_spectrum):
        self.grids = grids                  # GridManager: reference_grid() is the 1D LTE grid
        self.atom_provider = atom_provider  # () -> [atom paths] (the tool's own selector)
        self.add_spectrum = add_spectrum    # (spec, name) -> add the reference to the canvas

        saved = dict(store.get_section("abundcorr"))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        self.form1 = Form("corr1", [
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=PROFILES),
            Field("ew_min", "Min EW [mÅ]", kind="number", value=vals["ew_min"], min=0),
        ], on_change=self._on_field)
        self.form3 = Form("corr3", [
            Field("xaxis", "Scatter x-axis", kind="select", value=vals["xaxis"],
                  options=list(X_AXES)),
        ], on_change=self._on_field)

        super().__init__(ctx, setup=setup)
        # per-line error bars: the corrections carry a formal sigma(dA) since the
        # reference EW's fit-covariance error propagates through the CoG slope
        src = ColumnDataSource(dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[]))
        fig = self._build_fig(src)
        from bokeh.models import Whisker
        from canvas import ACCENT
        fig.add_layout(Whisker(base="x", upper="hi", lower="lo", source=src,
                               line_color=ACCENT, line_width=1.2))

    def refresh(self):
        pass

    # --- session persistence (whole computed results, saved into the session slot) ------
    def snapshot(self):
        """Serialise the computed corrections + selection for the session slot. Each result
        carries its own synthetic 1D LTE reference (``ref_spec``), so it is self-contained."""
        return {"results": [r.to_dict() for r in self.results],
                "selected": self.selected_uid}

    def restore(self, data):
        """Rebuild the corrections (and left-sidebar list) from a saved snapshot."""
        data = data or {}
        self.results = [ab.AbundanceResult.from_dict(d) for d in (data.get("results") or [])]
        uids = {r.uid for r in self.results}
        self.selected_uid = data.get("selected") if data.get("selected") in uids else None
        self.push_results()

    # --- state persistence -------------------------------------------------------------
    def _on_field(self, key, value):
        store.update_pref("abundcorr", key, value)
        if key == "xaxis":
            self.fig.xaxis.axis_label = X_AXES[str(value)][0]
            self._render_scatter()

    def _extra_action(self, kind, arg):
        if kind == "run-ref":
            self.run_reference()
        elif kind == "run-corr":
            self.run_corrections()

    # --- shared-machinery hooks -----------------------------------------------------------
    def _summarise(self, lines):
        return _summary(lines)

    @staticmethod
    def _summary_sigma(s):
        """The quotable uncertainty: the weighted formal error when the summary carries
        one, else the line-to-line scatter (legacy / unweighted results)."""
        sig = s.get("sigma_mean", float("nan"))
        if not np.isfinite(sig):
            sig = s.get("sigma_scatter", float("nan"))
        return sig

    def _result_label(self, res):
        # correction objects are named after the model atom they were computed from
        # (falls back to the persisted name for sessions restored without the atom path)
        return getattr(res, "atom_name", "") or res.name

    def _jump_decor(self, res, ln):
        # the best-fit NLTE synthesis at the inverted abundance, over the line window
        self.canvas.set_abund_overlay(ln.fit_spec, ln.window)
        if ln.fit_func is not None:
            self.canvas.set_abund_hatch(ln.fit_func, ln.window)
        # draw the fit window (inner, from auto_crop_core) and synth window (outer, HWHM-based)
        # bands at their true extents. Unlike AbundTool we do NOT snap the edges to the reference
        # spectrum's samples: the synthetic reference can be coarsely sampled, and snapping to a
        # coarse grid would misrepresent the auto_crop_core / HWHM extents the windows come from.
        self.canvas.set_abund_windows(ln.fit_window, ln.window)

    # --- helpers -----------------------------------------------------------------------
    @staticmethod
    def _species(ion):
        """(species, symbol) from a linelist ``ion`` label like ``"Fe I"``."""
        species = " ".join(str(ion).split())
        return species, (species.split()[0] if species else "")

    # --- step 1: reference + line selection --------------------------------------------
    def run_reference(self):
        ctx = self.setup.context()
        if ctx is None:
            self.notify(self.setup.reason or "Configure the target atmosphere (step 1).")
            return
        ref_grid = self.grids.reference_grid()
        if ref_grid is None:
            self.notify("Flag a 1D LTE reference grid in Preferences → Synthesis first.")
            return
        atoms = list(self.atom_provider())
        if not atoms:
            self.notify("Choose a model atom in the Corr panel (step 1) first.")
            return
        v = self.form1.values()
        profile, ew_min = str(v["profile"]), float(v["ew_min"])
        self._error = None
        self.notify("Computing the 1D LTE reference from the model atom…")
        loader = lambda token=None: self._compute_reference(
            ctx, ref_grid, atoms[0], profile, ew_min, token)
        self.ctx.submit("corr-reference", loader, self._on_reference_done)

    def _compute_reference(self, ctx, ref_grid, atom, profile, ew_min, token):
        try:
            prefs, target = ctx.prefs, ctx.target
            # the model atom alone defines the spectrum and its wavelength coverage, so no
            # window is imposed — the reference spans the atom's full synthesised range.
            ref_spec = synth.compute_reference(target, prefs, atom, ref_grid,
                                               cancel_token=token)
            llist = synth.read_reference_linelist(prefs)
            if llist is None or len(llist) == 0:
                self._error = "No linelist_out.txt was produced by the reference run."
                return None
            cand = llist.sort_values("lam0")
            if cand.empty:
                self._error = "The reference run produced no atom lines."
                return None
            lines = []
            for row in cand.itertuples():
                species, symbol = self._species(getattr(row, "ion", ""))
                line = {"lam0": float(row.lam0)}
                try:
                    lprof, win, fit_window, flags = ab.measure_line(
                        ref_spec, line, profile=profile)
                except Exception:
                    continue
                if "no_window" in flags or lprof is None or not np.isfinite(lprof.mAA):
                    continue
                ew, fwhm, func = lprof.mAA, float(lprof.fwhm), lprof.func
                # a fit region resolved from too few reference samples can't be fitted
                # reliably; flag it (and exclude it) so it doesn't skew the corrections.
                flags = list(flags)
                n_pts = len(ref_spec.crop(xmin=fit_window[0], xmax=fit_window[1]).values[0])
                if n_pts < MIN_FIT_POINTS:
                    flags.append("few_points")
                # keep the fit-covariance EW error: even on the noiseless synthetic
                # reference the residual-scaled covariance measures the profile/model
                # mismatch, and it propagates to a per-line sigma(dA) through the CoG
                # slope in step 2 — the formal error the weighted element mean uses.
                ew_err = float(lprof.width_err_mAA)
                # flag (and exclude) a line whose EW is more uncertain than the EW itself:
                # the fit-covariance error is the only EW uncertainty on the synthetic
                # reference, so ew_err > ew means the measurement is not trustworthy.
                if np.isfinite(ew_err) and ew > 0 and ew_err > ew:
                    flags.append("ew_unc")
                if not flags and ew < ew_min:
                    continue
                lines.append(ab.LineResult(
                    lam0=float(row.lam0), ion=species, e_low=float(getattr(row, "evi", np.nan)),
                    loggf=float(getattr(row, "loggf", np.nan)), ew=ew,
                    ew_err=ew_err, fwhm=fwhm,
                    mask=ab.LineMask(lam0=float(row.lam0), lmin=fit_window[0], lmax=fit_window[1],
                                     smin=win[0], smax=win[1], note=species),
                    fit_func=func, flags=tuple(flags), include=not flags))
            if not lines:
                self._error = f"No atom lines above {ew_min:.0f} mÅ in the reference."
                return None
            symbol = self._species(lines[0].ion)[1]
            res = ab.AbundanceResult(
                element=symbol, symbol=symbol, lines=lines, summary=_summary([]),
                stage="lines", profile=profile, uid=uuid.uuid4().hex[:8], relative=True)
            res.teff, res.logg, res.feh, res.vmic = target
            res.ref_spec = ref_spec
            res.atom_name = Path(atom).name
            res.name = res.atom_name
            return res
        except synth.CancelledError:
            raise
        except Exception as exc:
            logging.exception("Correction reference synthesis failed")
            self._error = str(exc)
            return None

    def _on_reference_done(self, result):
        if result is None:
            if self._error:
                self.notify(f"Reference failed: {self._error}")
            return
        self.results.append(result)
        self.selected_uid = result.uid
        # show the reference spectrum on the canvas so lines can be inspected/jumped to
        try:
            self.add_spectrum(result.ref_spec, f"{result.symbol} 1D LTE ref")
        except Exception:
            logging.exception("Could not add the reference spectrum to the canvas")
        self.set_step(1)
        self.select_result(result.uid)
        self.notify(f"{result.symbol}: {len(result.lines)} lines — review, then run corrections.")

    # --- step 2: per-line NLTE correction (EW match, streamed) -------------------------
    def run_corrections(self):
        res = self.selected()
        if res is None:
            self.notify("Select a correction object first.")
            return
        ctx = self.setup.context()
        if ctx is None:
            self.notify(self.setup.reason or "Configure the target atmosphere (step 1).")
            return
        atoms = list(self.atom_provider())
        if not atoms:
            self.notify("Choose a model atom in the Corr panel (step 1).")
            return
        self._error = None
        self.notify(f"Computing NLTE corrections for {res.symbol}…")
        loader = lambda token=None: self._corr_lines(res, ctx, atoms[0], token)
        self.ctx.submit(f"corr-synth:{res.symbol}", loader, self._on_corr_done)

    def _corr_lines(self, res, ctx, atom, token):
        try:
            prefs, grid = ctx.prefs, ctx.grid
            target, atmos_file = ctx.target, ctx.atmos_file
            symbol = res.symbol
            ladder = np.linspace(-CORR_SPAN, CORR_SPAN, CORR_N)

            def _rv(a):
                return {"mode": "NLTE", "only_atom": True, "atoms": [atom],
                        "abundances": f"{symbol}={a:.4f}", "relative_abundances": True}

            # set_spectrum=False leaves spectrum_params without R/daa, so only the aa_blue/aa_red
            # window is passed: the Fortran samples the line on the model atom's own frequency
            # grid (matching the 1D LTE reference) and crops it to the synth window, rather than
            # resampling onto a uniform R grid.
            def nlte_synth(a, window):
                return synth.compute(target, prefs, _rv(a), window[0] - SYNTH_PAD_AA,
                                     window[1] + SYNTH_PAD_AA, grid=grid, atmos_file=atmos_file,
                                     set_spectrum=False, cancel_token=token)

            # the whole abundance ladder for one line synthesised in a single call, so every
            # point's dispatch run is fanned across one shared process pool instead of being
            # issued one after another (results are identical to the serial nlte_synth loop).
            def batch_nlte_synth(a_list, window):
                return synth.compute_batch(target, prefs, [_rv(a) for a in a_list],
                                           window[0] - SYNTH_PAD_AA, window[1] + SYNTH_PAD_AA,
                                           grid=grid, atmos_file=atmos_file,
                                           set_spectrum=False, cancel_token=token)

            for ln in res.lines:
                if token is not None and token.cancelled:
                    raise synth.CancelledError()
                if not ln.include:
                    continue
                try:
                    a_grid, ew_grid = ab.build_cog(
                        None, ladder, window=ln.fit_window,
                        profile=res.profile, x0=ln.lam0,
                        batch_synth_fn=lambda a_list, _w=ln.window: batch_nlte_synth(a_list, _w))
                    # propagate the reference EW's fit-covariance error through the local
                    # CoG slope: sigma(dA) = sigma_EW * |dA/dEW|. All of it is fit/model
                    # error (no continuum term on a synthetic reference).
                    ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                    inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
                    ln.a = float(inv["a"])
                    ln.sigma_a = float(inv["sigma_a"])
                    ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
                    ln.a_grid, ln.ew_grid = tuple(a_grid), tuple(ew_grid)
                    if inv["out_of_range"]:
                        ln.flags = ("out_of_range",)
                        ln.include = False
                    else:
                        # overlay the best-fit NLTE line — synthesised at the inverted
                        # abundance — so jumping to a row shows its fit against the
                        # reference, mirroring the abundance tool. A failed overlay must
                        # not discard the (successful) correction, so guard it separately.
                        try:
                            ln.fit_spec = nlte_synth(ln.a, ln.window)
                        except synth.CancelledError:
                            raise
                        except Exception:
                            logging.exception("Correction: overlay synth failed %.3f", ln.lam0)
                except synth.CancelledError:
                    raise
                except Exception:
                    logging.exception("Correction: line %.3f failed", ln.lam0)
                    ln.flags = ("synth_failed",)
                self.ctx.schedule(lambda r=res: self._stream_row(r))
            return res
        except synth.CancelledError:
            raise
        except Exception as exc:
            logging.exception("Correction synthesis failed")
            self._error = str(exc)
            return None

    def _on_corr_done(self, result):
        if result is None:
            if self._error:
                self.notify(f"Corrections failed: {self._error}")
            return
        result.stage = "corr"
        result.summary = _summary(result.lines)
        s = result.summary
        result.name = self._result_label(result)
        self.push_results()
        if self.selected() is result:
            self._render_selected()
        if np.isfinite(s["a_mean"]):
            self.notify(f"{result.symbol}: ⟨ΔA⟩ = {s['a_mean']:+.3f} ± "
                        f"{self._summary_sigma(s):.3f} (n={s['n']}) — see the scatter in step 3.")

    # --- rendering ----------------------------------------------------------------------
    def _render_scatter(self):
        res = self.selected()
        inc = [ln for ln in res.lines if ln.include and np.isfinite(ln.a)] if res else []
        if not inc:
            self.src.data = dict(x=[], a=[], lo=[], hi=[], lam=[], ion=[])
            self.mean_line.visible = self.mean_band.visible = False
            return
        xf = X_AXES[str(self.form3.value("xaxis"))][1]
        xs = [xf(ln) for ln in inc]
        # zero-length whiskers for lines without a formal error (legacy sessions)
        sig = [ln.sigma_a if np.isfinite(ln.sigma_a) else 0.0 for ln in inc]
        los = [ln.a - s for ln, s in zip(inc, sig)]
        his = [ln.a + s for ln, s in zip(inc, sig)]
        self.src.data = dict(x=xs, a=[ln.a for ln in inc], lo=los, hi=his,
                             lam=[ln.lam0 for ln in inc], ion=[ln.ion for ln in inc])
        self._crop_scatter(xs, los + his)
        s = res.summary
        if np.isfinite(s.get("a_mean", np.nan)):
            self.mean_line.location = s["a_mean"]
            self.mean_line.visible = True
            band = self._summary_sigma(s)
            band = band if np.isfinite(band) else 0.0
            self.mean_band.top = s["a_mean"] + band
            self.mean_band.bottom = s["a_mean"] - band
            self.mean_band.visible = True
        else:
            self.mean_line.visible = self.mean_band.visible = False

    def results_html(self):
        if not self.results:
            return ""
        rows = self._results_rows("Δ", " corr-result-item", "data-corr-result",
                                  "data-corr-remove")
        return ('<div class="section-header result-section-head">Corrections</div>'
                + "".join(rows))

    def table_html(self, res):
        head = self._table_head("corr-mask-chk", (
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), ("ΔA", "a"), ("flags", None)))
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
            if self.show_toggled_only and not ln.include:
                continue
            chk = " checked" if ln.include else ""
            cls = "" if ln.include else " line-excluded"
            ew = self._fmt(ln.ew, ln.ew_err, digits=1)
            da = self._fmt(ln.a, ln.sigma_a if res.stage == "corr" else None, plus=True)
            fl = ", ".join(ln.flags)
            body.append(
                f'<tr class="corr-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="corr-line-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'<td>{da}</td>'
                f'<td class="abund-flags">{html.escape(fl)}</td></tr>')
        foot = ""
        if res.stage == "corr" and np.isfinite(res.summary.get("a_mean", np.nan)):
            s = res.summary
            if np.isfinite(s.get("sigma_mean", np.nan)):
                budget = (f'formal {s["sigma_mean"]:.3f} · '
                          f'scatter {s["sigma_scatter"]:.3f} · n={s["n"]}')
            else:
                budget = f'line-to-line σ · n={s["n"]}'
            foot = (f'<div class="abund-summary">⟨ΔA({html.escape(res.symbol)})⟩ = '
                    f'{s["a_mean"]:+.3f} ± {self._summary_sigma(s):.3f}'
                    f'<span class="abund-budget"> ({budget})</span></div>')
        return (f'<table class="linelist-table abund-table corr-lines">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body)}</tbody></table>{foot}')
