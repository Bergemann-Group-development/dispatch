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

All three steps are :mod:`tumult.abundance` calls — :func:`~tumult.abundance.measure_lines`,
:func:`~tumult.abundance.run_abundance`, :func:`~tumult.abundance.finalize_result` — driven
here on the setup's atmosphere, with a ``progress`` callback streaming each line into the
table as it lands. What is left in this module is the wizard's own machinery, shared with the
corrections tool through :mod:`abundwizard`, and the rendering.
"""
from __future__ import annotations

import html
import logging
import astropy.units as u
import numpy as np
from bokeh.models import ColumnDataSource, Whisker

import abundance as ab
import store
import synth
from abundwizard import PROFILES, X_AXES, AbundanceWizard
from canvas import ACCENT
from controls import Field, Form

DEFAULTS = {
    "profile": "gauss", "snr": 100.0, "cont_pct": 0.5, "xaxis": "evi",
}


class AbundanceTool(AbundanceWizard):
    name = "abundance"
    prefix = "abund"
    kind = "result"
    btn_selector = "#abundance-btn"
    empty_table_msg = "No abundance result selected"
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
        # step 1: what to measure (independent of stellar parameters). The EW error inputs
        # (SNR / continuum-%) live here too, so each line's EW uncertainty is known at
        # selection; a line whose photon EW error exceeds its EW is flagged automatically.
        # The lines to measure come from a line-mask set (Objects → Open File → New Line
        # Mask), so the element and the candidate search live with the generator that made
        # it; step 1 only measures the windows it is given, on the chosen reference.
        self.form1 = Form("abund1", [
            Field("maskset", "Line mask", kind="select", value="", options=[]),
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=PROFILES),
            Field("snr", "SNR (observation)", kind="number", value=vals["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=vals["cont_pct"], min=0, step=0.1),
        ], on_change=self._on_field)
        self.maskset_uid = saved.get("maskset_uid")
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
        self._sync_masksets()                       # step-1 line-mask selector

    def refresh(self):
        """Keep the step-1 reference + line-mask selectors in sync as objects come and go."""
        if self.active:
            self._ref_bridge.value = self.ref_html()
            self._sync_masksets()

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
            ln.ew, ln.ew_err = prof.maa, prof.width_err_maa
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
        if key == "maskset":
            self.maskset_uid = value or None
            store.update_pref("abundance", "maskset_uid", self.maskset_uid or "")
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

    # --- step 1: diagnostic-line selection ---------------------------------------------
    def run_lines(self):
        ref_uid, ref_spec, err = self._ref()
        if err:
            self.notify(err)
            return
        mset = self.selected_maskset()
        if mset is None:
            self.notify("Make a line mask first — Objects → Open File → New Line Mask.")
            return
        if not mset.element:
            self.notify(f"Set the element of “{mset.name}” before measuring it.")
            return
        if not any(m.include for m in mset.masks):
            self.notify(f"Every window in “{mset.name}” is switched off.")
            return
        v = self.form1.values()
        self._error = None
        self.notify(f"Measuring {mset.element} windows on the reference…")
        loader = lambda token=None: self._measure_lines(
            ref_uid, ref_spec, mset, str(v["profile"]), float(v["snr"]) or None,
            float(v["cont_pct"]))
        self.ctx.submit(f"abund-lines:{mset.element}", loader, self._on_lines_done)

    def _measure_lines(self, ref_uid, ref_spec, mset, profile, snr, cont_pct):
        """Measure the set's switched-on windows (:func:`tumult.abundance.measure_lines`)."""
        try:
            res = ab.measure_lines(mset, ref_spec, linelist=self.master_path,
                                   profile=profile, snr=snr, cont_pct=cont_pct)
            res.ref_uid = ref_uid
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
        """Fit every included line (:func:`tumult.abundance.run_abundance`), streaming rows."""
        obs = self.get_observation(getattr(res, "ref_uid", None))
        if obs is None:
            self._error = "The reference observation for this result is no longer loaded."
            return None
        try:
            return ab.run_abundance(
                ctx, res, obs, cancel_token=token,
                progress=lambda r: self.ctx.schedule(lambda r=r: self._stream_row(r)))
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
        ab.finalize_result(res, spec, float(v["snr"]) or None, float(v["cont_pct"]))
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
        head = self._table_head("abund-mask-chk", (
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), ("blend", "blend"),
            ("A(X)", "a"), ("vshift", "vshift"), ("vdop", "vdop"), ("flags", None)))
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
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
        return (f'<table class="linelist-table abund-table abund-lines">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body)}</tbody></table>{foot}')
