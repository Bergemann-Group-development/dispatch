"""The Abundance-Correction tool: a 3-step wizard for **NLTE abundance corrections vs 1D LTE**.

Toolbar **Corr** opens a right-sidebar panel whose controls follow a step switcher shown in a
footer inside the centre panel. It reuses the abundance pipeline numerics
(:mod:`abundance` / :mod:`tumult.abundance`), the :class:`synthsetup.SynthSetup` config cluster,
the machinery shared with the abundance wizard (:mod:`abundwizard`), and the app-global grid
registry — one grid of which is flagged as the **1D LTE reference**.

1. **Reference & lines.** Pick stellar parameters (or a direct atmosphere file — the NLTE target)
   and a single model atom. "Compute reference" runs
   :func:`tumult.abundance.measure_reference_lines`, which interpolates the flagged 1D LTE grid
   at those parameters, synthesises the reference **straight from the model atom** (no line
   list), and measures every candidate line from that run's ``linelist_out.txt``.
2. **NLTE corrections.** :func:`tumult.abundance.run_nlte_corrections` builds an NLTE curve of
   growth per line on the *picked* atmosphere (model atom only) and inverts it to the abundance
   reproducing the reference EW. The reference sits at ``[X/Fe]=0``, so that abundance **is**
   the correction ΔA = A(NLTE) − A(LTE). Rows stream into the table via the ``progress`` hook.
3. **Scatter.** Plot the line-to-line scatter of the corrections (mean ± population σ).

Both stages are headless: what is left here is the form handling, submission and rendering.
"""
from __future__ import annotations

import html
import logging

import numpy as np
from bokeh.models import ColumnDataSource

import abundance as ab
import store
import synth
from abundwizard import PROFILES, X_AXES, AbundanceWizard
from controls import Field, Form
# The workflow itself is headless (:mod:`tumult.abundance`); what is left here is the
# wizard: collecting the form values, submitting the two stages, and rendering their results.
from tumult.abundance import (CorrectionOpts, correction_summary,
                              measure_reference_lines, run_nlte_corrections)

DEFAULTS = {"profile": "gauss", "ew_min": 5.0, "xaxis": "evi"}


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
        return correction_summary(lines)

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
        """Stage 1, on a worker thread: the headless call plus this tool's error reporting."""
        try:
            return measure_reference_lines(
                ctx, atom, ref_grid=ref_grid,
                opts=CorrectionOpts(profile=profile, ew_min=ew_min),
                cancel_token=token)
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
        """Stage 2, on a worker thread: rows stream into the table as each line lands."""
        try:
            return run_nlte_corrections(
                ctx, res, atom, opts=CorrectionOpts(profile=res.profile),
                progress=lambda r: self.ctx.schedule(lambda _r=r: self._stream_row(_r)),
                cancel_token=token)
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
