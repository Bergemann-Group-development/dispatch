"""The Abundance-Correction tool: **NLTE abundance corrections vs 1D LTE**, line by line.

The entry point is a model-atom object (:mod:`atomtool`): its panel's **NLTE corrections…**
action opens this tool's modal with the atom fixed. The modal holds the NLTE target
atmosphere and stellar parameters (a :class:`synthsetup.SynthSetup` cluster), the measurement
settings, and states which registered grid is flagged as the **1D LTE reference** (set in
Preferences → Synthesis). Its one **Compute** button runs both stages as a chained job:

1. **Reference.** :func:`tumult.abundance.measure_reference_lines` interpolates the flagged
   1D LTE grid at the target parameters, synthesises the reference **straight from the model
   atom** (no line list) and measures every candidate line from that run's
   ``linelist_out.txt``. The reference spectrum lands on the canvas and the result object
   appears under its atom in the left sidebar as soon as this stage finishes.
2. **Corrections.** :func:`tumult.abundance.run_nlte_corrections` builds an NLTE curve of
   growth per line on the picked atmosphere (model atom only) and inverts it to the abundance
   reproducing the reference EW. The reference sits at ``[X/Fe]=0``, so that abundance **is**
   the correction ΔA = A(NLTE) − A(LTE). Rows stream into the table via the ``progress`` hook.

The result view has two steps: **Lines** (the reference spectrum over the per-line table) and
**Scatter** (the line-to-line scatter of the corrections, mean ± population σ). Both stages
are headless: what is left here is the modal, the submission and the rendering.
"""
from __future__ import annotations

from copy import deepcopy

import html
import logging
from pathlib import Path

import numpy as np
from bokeh.models import ColumnDataSource

from tumult import abundance as ab
import store
from abundwizard import PROFILES, X_AXES, AbundanceWizard
from controls import Field, Form
from htmlbridge import pulse
from tumult.abundance import (CorrectionOpts, correction_summary,
                              measure_reference_lines, run_nlte_corrections)

DEFAULTS = {"profile": "gauss", "ew_min": 5.0, "xaxis": "evi"}


class AbundCorrTool(AbundanceWizard):
    name = "abundcorr"
    prefix = "corr"
    kind = "corr"
    steps = 2                        # lines, scatter
    empty_table_msg = "No correction selected"
    y_axis_label = "ΔA (NLTE − LTE)"
    final_stage = "corr"

    def __init__(self, ctx, *, setup, grids, add_spectrum):
        self.grids = grids                  # GridManager: reference_grid() is the 1D LTE grid
        self.add_spectrum = add_spectrum    # (spec, name) -> add the reference to the canvas
        self.on_results = lambda: None      # installed by main: repaints the atom list

        saved = dict(store.get_section("abundcorr"))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # the measurement settings live in the modal; the scatter axis in the result panel
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
        self._modal = ctx.bridges.dialog("corr-modal", name="corr-modal-bridge")
        self._context_bridge = ctx.bridges.html("corr-context")
        self._pending_atom = None           # the AtomObject the open modal computes for
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

    # --- results by atom (the left-sidebar nesting) ------------------------------------
    def results_for(self, atom_path):
        """The results computed from ``atom_path`` (``None``: those whose atom is not
        loaded, i.e. whose path matches no listed atom)."""
        if atom_path is None:
            known = set(self._known_atom_paths())
            return [r for r in self.results if r.atom_path not in known]
        return [r for r in self.results if r.atom_path == str(atom_path)]

    def _known_atom_paths(self):
        """Installed by main: the paths of the loaded atom objects."""
        return []

    def result_rows(self, results, nested=False):
        return self._results_rows("Δ", " corr-result-item", "data-corr-result",
                                  "data-corr-remove", results=results, nested=nested)

    def push_results(self):
        self.on_results()

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
        if kind == "run":
            self.run()

    # --- the modal ----------------------------------------------------------------------
    def open_for(self, atom_obj):
        """Open the modal to compute corrections from ``atom_obj`` (the atom panel action)."""
        self._pending_atom = atom_obj
        if not store.get_section(self.setup.section).get("_seeded"):
            self.setup.inherit_from()
            store.update_pref(self.setup.section, "_seeded", True)
        self.setup.render()
        self.form1.render()
        self._push_modal_context()
        pulse(self._modal, "open")

    def _push_modal_context(self):
        self._context_bridge.value = self.context_html()

    def context_html(self):
        """The modal's header block: the fixed atom and the flagged 1D LTE reference."""
        obj = self._pending_atom
        if obj is None:
            return '<div class="hint">No model atom selected.</div>'
        nk, nl, _nc = obj.counts()
        ref = self.grids.reference_grid()
        if ref is None:
            ref_html = ('<div class="hint hint-warn">No 1D LTE reference grid — register a '
                        'grid and flag it in Preferences → Synthesis.</div>')
        else:
            ref_name = html.escape(ref.get("name") or ref["id"])
            ref_html = (f'<div class="hint">1D LTE reference: <b>{ref_name}</b>, interpolated '
                        'at the parameters below (flagged in Preferences → Synthesis).</div>')
        return (f'<div class="hint">NLTE corrections for <b>{html.escape(obj.symbol)}</b> '
                f'({html.escape(obj.name)}) — {nl} bound-bound transitions between {nk} '
                'levels. The atom defines the spectrum and its wavelength range.</div>'
                + ref_html)

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
        # (falls back to the persisted name for sessions saved without the atom path)
        return Path(res.atom_path).name if res.atom_path else res.name

    def _jump_decor(self, res, ln):
        # the best-fit NLTE synthesis at the inverted abundance, over the line window
        self.canvas.set_abund_overlay(ln.fit_spec, ln.window)
        if ln.fit_func is not None:
            self.canvas.set_abund_hatch(ln.fit_func, ln.window)
        # highlight the fit window (inner, from auto_crop_core) and synth window (outer,
        # HWHM-based) at their true extents. Unlike AbundTool we do NOT snap the edges to the
        # reference spectrum's samples: the synthetic reference can be coarsely sampled, and
        # snapping to a coarse grid would misrepresent the auto_crop_core / HWHM extents the
        # windows come from.
        self.canvas.highlight_windows([ln.mask], active=0)

    # --- the modal's run: reference, then corrections, as one chained job ---------------
    def run(self):
        obj = self._pending_atom
        if obj is None:
            return
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason or "Configure the target atmosphere.")
            return
        ref_grid = self.grids.reference_grid()
        if ref_grid is None:
            self.notify("Flag a 1D LTE reference grid in Preferences → Synthesis first.")
            return
        pulse(self._modal, "close")
        v = self.form1.values()
        profile, ew_min = str(v["profile"]), float(v["ew_min"])
        self.notify(f"Computing the 1D LTE reference for {obj.symbol}…")
        ref_grid = deepcopy(ref_grid)
        atom = obj.path
        loader = lambda _res, token, progress: measure_reference_lines(
            request, atom, ref_grid=ref_grid,
            opts=CorrectionOpts(profile=profile, ew_min=ew_min), cancel_token=token)
        self.jobs.submit(f"corr-reference:{obj.symbol}", loader,
                         lambda result: self._on_reference_done(result, request))

    def _on_reference_done(self, result, request):
        self.results.append(result)
        self.selected_uid = result.uid
        # show the reference spectrum on the canvas so lines can be inspected/jumped to
        try:
            self.add_spectrum(result.ref_spec, f"{result.symbol} 1D LTE ref")
        except Exception:
            logging.exception("Could not add the reference spectrum to the canvas")
        self.set_step(1)
        self.select_result(result.uid)
        self.notify(f"{result.symbol}: {len(result.lines)} reference lines — computing "
                    "NLTE corrections…")
        self._run_corrections(result, request)

    def _run_corrections(self, res, request):
        atom = res.atom_path
        loader = lambda working, token, progress: run_nlte_corrections(
            request, working, atom, opts=CorrectionOpts(profile=working.profile),
            progress=progress, cancel_token=token)
        self.jobs.submit(f"corr-synth:{res.symbol}", loader, self._on_corr_done, result=res)

    def _on_corr_done(self, result):
        s = result.summary
        result.name = self._result_label(result)
        self.push_results()
        if self.selected() is result:
            self._render_selected()
        if np.isfinite(s["a_mean"]):
            self.notify(f"{result.symbol}: ⟨ΔA⟩ = {s['a_mean']:+.3f} ± "
                        f"{self._summary_sigma(s):.3f} (n={s['n']}) — see the Scatter step.")

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
