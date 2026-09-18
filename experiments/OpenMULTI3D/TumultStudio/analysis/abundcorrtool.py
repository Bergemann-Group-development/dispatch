"""The Abundance-Correction tool: **NLTE abundance corrections vs 1D LTE**, line by line.

The entry point is a model-atom object (:mod:`atomtool`): its panel's **NLTE corrections…**
action opens this tool's modal with the atom fixed. The modal holds the NLTE target
atmosphere and stellar parameters (a :class:`synthsetup.SynthSetup` cluster), the measurement
settings, and states which registered grid is flagged as the **1D LTE reference** (set in
Preferences → Atmosphere grids). Its one **Compute** button runs both stages as one durable
``abundance.nlte_corrections`` workflow (:mod:`tumult.synthesis.tasks`, through
``ResultJobs.submit_operation``) on the setup's request with the atom as requested
departures and the flagged grid as ``reference_atmosphere``. It outlives the tab: the
stage-1 result arrives as the record's first progress snapshot, every correction round as
the next ones.

1. **Reference.** :func:`tumult.abundance.measure_reference_lines` interpolates the flagged
   1D LTE grid at the target parameters and composition, synthesises the reference
   **straight from the model atom** (no line list) and measures every candidate line from that run's
   ``linelist_out.txt``. By default (*Line grid from opacity profiles*) the bound-bound
   wavelength grid and the per-line windows are generated first from the atom's line
   opacities on the reference model nearest the target (:mod:`tumult.line_quadrature`);
   those windows are the line masks, and stage 2 runs on the same grid. Without the line
   list, a window that overlaps a neighbouring atom line's is narrowed to the reference
   flux maxima nearest the centre, so a blended line is integrated over its own profile.
   The result object appears under its atom in the left sidebar as soon as this stage
   finishes; the reference spectrum itself is drawn only while this result's correction view
   is shown (:meth:`_render_reference`), not as a sidebar spectrum of its own.
   With **Use line list** on, the configured synthesis line list is synthesised
   with the atom (reference and every NLTE rung), and the reference run's ``linelist_out.txt``
   is what the canvas annotates while the result is shown; off, the canvas shows no line
   annotations for it, since the spectrum holds nothing but the atom's lines.
2. **Corrections.** :func:`tumult.abundance.run_nlte_corrections` builds one NLTE curve of
   growth for all lines on the picked atmosphere: each rung is a full NLTE solve at an
   abundance offset on the preferences' ``abund_step`` lattice within the modal's **CoG
   limit**, synthesising the atom's whole range, and the rungs are added round by round
   until every line's inverted correction is bracketed by two computed rungs. The rungs
   are offsets from the element's abundance in the setup's composition, the reference's
   own, so the inverted offset **is** the correction ΔA = A(NLTE) − A(LTE). Rows stream
   into the table via the ``progress`` hook after every round.

The result view has two steps: **Lines** (the reference spectrum, with the selected line's
curve-of-growth inversion beside it, over the per-line table) and **Scatter** (the
line-to-line scatter of the corrections against excitation energy and reduced EW, mean ± σ).
On the Lines step the jumped-to line's fit window is the accented box and its edges can be
dragged, within the synth window the atom fixes; :func:`tumult.abundance.refit_line_window`
re-reads the reference EW and the stored rungs over the new window, so the correction
follows the edit without another solve. An edited line's row offers **Reset window**
(:func:`tumult.abundance.reset_line_window`), which returns to the window the line was
computed over by the same re-read. Both stages are headless: what is left here is the
modal, the submission and the rendering.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

import numpy as np
from bokeh.models import ColumnDataSource

from tumult import abundance as ab
from ..workspace import naming
from ..workspace import store
from .diagnostics import empty_points
from .resultview import ResultView, format_value
from ..preferences.synthprefs import atmosphere_grid
from tumult.departures import DepartureRequest
from ..spectrum.canvas import ACCENT
from .cogpanel import CogPanel
from ..ui.controls import Field, Form
from ..spectrum.ew import FIT_PROFILES, profile_options
from ..ui.htmlbridge import pulse
from tumult.abundance import (CorrectionOpts, correction_summary, reduced_ew,
                              refit_line_window, reset_line_window, window_edited)
from ..ui.units import to_range

import astropy.units as u

# the trapezoidal width integrates each line over its atom window, so it needs the line
# grid; the fitted profiles work either way (see tumult.abundance.CorrectionOpts)
CORR_PROFILES = ["trapz"] + list(FIT_PROFILES)
DEFAULTS = {"profile": "trapz", "ew_min": 5.0, "limit": 1.0, "line_grid": True,
            "linelist": False}
REF_COLOR = "#7f8a99"        # the reference spectrum's own line: muted, so the ACCENT-drawn
                              # fit overlays and windows still read as what the tool is about


class AbundCorrTool(ResultView):
    name = "abundcorr"
    prefix = "corr"
    kind = "corr"
    steps = 2                        # lines, scatter
    empty_table_msg = "No correction selected"
    y_axis_label = "ΔA (NLTE − LTE)"
    final_stage = "corr"
    extra_sort_keys = {"rew": lambda ln: reduced_ew(ln.ew, ln.lam0)}
    extra_sort_asc = {"rew": False}

    def __init__(self, ctx, *, setup, grids):
        self.grids = grids                  # GridManager: reference_grid() is the 1D LTE grid
        self.on_results = lambda: None      # installed by main: repaints the atom list

        saved = dict(store.get_section("abundcorr"))
        vals = {**DEFAULTS, **{k: saved[k] for k in DEFAULTS if k in saved}}
        # the measurement settings live in the modal
        self.form1 = Form("corr1", [
            Field("profile", "EW profile", kind="select", value=vals["profile"],
                  options=profile_options(CORR_PROFILES)),
            Field("ew_min", "Min EW [mÅ]", kind="number", value=vals["ew_min"], min=0),
            Field("limit", "Curve-of-growth half-span [dex]", kind="number", value=vals["limit"], min=0,
                  step=0.25),
            Field("line_grid", "Line grid from opacity profiles", kind="checkbox",
                  value=vals["line_grid"]),
            Field("linelist", "Use line list", kind="checkbox",
                  value=vals["linelist"]),
        ], on_change=self._on_field)

        super().__init__(ctx, setup=setup)
        self._modal = ctx.bridges.dialog("corr-modal", name="corr-modal-bridge")
        self._context_bridge = ctx.bridges.html("corr-context")
        self._pending_atom = None           # the AtomObject the open modal computes for
        # The reference spectrum: a renderer this view owns, not a canvas object of its own
        # (it is not user data, has no fit or window of its own, and should not persist in
        # the sidebar once the correction result it belongs to is not being looked at). It
        # does count towards the main figure's y auto-fit, unlike a tool's decorative
        # overlays - it *is* the spectrum being viewed while this result has focus.
        self.ref_source = ColumnDataSource(dict(x=[], y=[]))
        self.ref_line = self.canvas.main.line(
            "x", "y", source=self.ref_source, color=REF_COLOR, line_width=1.2,
            line_alpha=0.8, visible=False)
        self.canvas.include_in_y_fit(self.ref_line)
        # per-line error bars: the corrections carry a formal sigma(dA) since the
        # reference EW's fit-covariance error propagates through the CoG slope
        self.build_panels()
        self.panels.add_whiskers(self.panels.src, ACCENT)
        self.cog = CogPanel("corr_cog", crop=self.panels.crop, fmt=format_value)

    def _clear_line_decor(self):
        super()._clear_line_decor()
        self.cog.clear()

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
        return self._results_rows("Δ", results=results, nested=nested)

    def push_results(self):
        self.on_results()

    # --- session persistence ------------------------------------------------------------
    # the envelope is ResultView's, the factory and the relink :mod:`results`'. A correction
    # carries its own reference spectrum (``ref_spec``) and names its atom by path; a result
    # whose atom is not restored stays listed in the trailing Corrections section (see
    # results_for).

    # --- state persistence -------------------------------------------------------------
    def _on_field(self, key, value):
        pass  # The calculation draft is persisted on submission.

    def _extra_action(self, kind, arg):
        if kind == "run":
            self.run()
        elif kind == "reset-window":
            self.reset_window(int(arg))

    # --- the modal ----------------------------------------------------------------------
    def open_for(self, atom_obj):
        """Open the modal to compute corrections from ``atom_obj`` (the atom panel action)."""
        self._pending_atom = atom_obj
        self.draft.begin()
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
                        'grid and flag it in Preferences → Atmosphere grids.</div>')
        else:
            ref_name = html.escape(ref.get("name") or ref["id"])
            ref_html = (f'<div class="hint">1D LTE reference: <b>{ref_name}</b>, interpolated '
                        'at the parameters below (flagged in Preferences → Atmosphere grids).</div>')
        return (f'<div class="hint">NLTE corrections for <b>{html.escape(obj.symbol)}</b> '
                f'({html.escape(obj.name)}) — {nl} bound-bound transitions between {nk} '
                'levels. The atom defines the spectrum and its wavelength range. With the '
                'line grid on, its bound-bound sampling and each line\'s window are built '
                'from the line opacities on the reference model nearest the target, and '
                'those windows are the line masks; without the line list, overlapping '
                'windows are cut at the reference flux maxima around each line.</div>'
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

    def _result_summary(self, res):
        s = res.summary
        if res.stage == "corr" and np.isfinite(s.get("a_mean", np.nan)):
            return f"⟨ΔA⟩ = {s['a_mean']:+.3f} ± {self._summary_sigma(s):.3f}"
        return f"{len(res.lines)} reference lines"

    def _jump_decor(self, res, ln):
        self.jumped = next((i for i, l in enumerate(res.lines) if l is ln), None)
        self._render_line_decor(ln)

    def _render_line_decor(self, ln):
        """The jumped-to line's overlays: the rungs interpolated at its correction, the
        EW-area hatch of the reference fit, its windows and its curve of growth. The fit
        window is the accented box and can be dragged by its edges (:meth:`_on_fit_drag`);
        the synth window is the atom's own and only displayed. The edges are shown at their
        true extents, not snapped to the reference's samples, which can be coarse."""
        if ln.fit_spec is not None:
            self.canvas.set_abund_overlay(ln.fit_spec, ln.window)
        else:
            self.canvas.clear_abund_overlay()
        if ln.fit_func is not None:
            self.canvas.set_abund_hatch(ln.fit_func, ln.fit_window)
        else:
            self.canvas.clear_abund_hatch()
        self.canvas.highlight_windows([ln.mask], active=0, editable=True)
        self.cog.render(ln)

    def _on_fit_drag(self, x0, x1):
        """An edge of the accented fit window was dragged (display units): re-read the
        line over the new window, clamped to its synth window, and repaint everything that
        shows it - the table row, the scatter, the CoG panel and the mean."""
        res = self.selected()
        if res is None or self.jumped is None or not self.is_focused():
            return
        try:
            ln = res.lines[self.jumped]
        except IndexError:
            return
        lo, hi = to_range(x0, x1, self.canvas.xunit, u.AA)
        refit_line_window(res, ln, min(lo, hi), max(lo, hi))
        self._line_window_changed(res, ln)

    def reset_window(self, index):
        """**Reset window**: return line ``index`` to the fit window it was computed over."""
        res = self.selected()
        if res is None or not 0 <= index < len(res.lines):
            return
        ln = res.lines[index]
        if not window_edited(ln):
            return
        reset_line_window(res, ln)
        self.jumped = index
        self._line_window_changed(res, ln)

    def _line_window_changed(self, res, ln):
        """Repaint everything that shows a line's window: its overlays, the table row,
        the scatter, the CoG panel and the mean."""
        self._render_line_decor(ln)
        self.push_table(self.table_html(res))
        self._render_scatter()
        self.push_results()

    # --- line annotations: the reference run's own dump, or none ------------------------
    def select_result(self, uid):
        super().select_result(uid)
        self.jumped = None
        self.canvas.on_window_resize = self._on_fit_drag
        self._attach_annotations()
        self._render_reference()

    def blur(self):
        super().blur()
        self.canvas.on_window_resize = None
        self._detach_annotations()
        self._hide_reference()

    def deactivate(self):
        super().deactivate()
        self.canvas.on_window_resize = None
        self._detach_annotations()
        self._hide_reference()

    def remove_result(self, uid):
        placed = super().remove_result(uid)   # re-selects (and re-attaches) when one is left
        if self.selected() is None:
            self.canvas.on_window_resize = None
            self._detach_annotations()
            self._hide_reference()
        return placed

    # --- the reference spectrum: a renderer owned by this view -----------------------------
    def _render_reference(self):
        """Draw the selected result's 1D LTE reference on the main figure, or hide it.

        The reference is not user data and carries no fit or window of its own, so it is not
        a canvas object (no sidebar entry, nothing to toggle or remove independently) - it is
        only ever a backdrop for whichever correction result owns the centre view."""
        res = self.selected()
        ref_spec = getattr(res, "ref_spec", None) if res is not None else None
        if not (self.active and self.is_focused()) or ref_spec is None:
            return self._hide_reference()
        spec = ref_spec.decimate(self.canvas.n_main).convert_units(
            xunit=self.canvas.xunit, yunit=self.canvas.yunit)
        x, y = spec.values
        self.ref_source.data = dict(x=np.asarray(x, dtype=float), y=np.asarray(y, dtype=float))
        self.ref_line.visible = True

    def _hide_reference(self):
        self.ref_line.visible = False

    def _attach_annotations(self):
        """Own the line-position overlay while a correction is shown: the transitions the
        reference run actually used (its ``linelist_out.txt``, every species with HOP's
        strengths) when the line list was included, and nothing at all otherwise — an
        atom-only reference has no other species to name, and the default whole-list overlay
        would suggest lines the spectrum does not contain."""
        annot = self.ctx.annot
        res = self.selected()
        if annot is None or res is None:
            return
        lines = self._dump_lines(res)
        annot.attach(self.name, lines=lines, title=res.symbol if len(lines) else "")

    def _detach_annotations(self):
        if self.ctx.annot is not None:
            self.ctx.annot.detach(self.name)

    def _dump_lines(self, res):
        """``lam0``/``ion``/``strength`` of the result's reference dump, cached per result;
        an empty frame when the run had no line list or the dump is gone."""
        import pandas as pd
        cache = self.__dict__.setdefault("_dump_cache", {})
        dump = str(getattr(res, "linelist_dump", "") or "")
        if dump not in cache:
            frame = pd.DataFrame({"lam0": [], "ion": [], "strength": []})
            if dump and Path(dump).is_file():
                try:
                    from tumult.linelist import LineList
                    llist = LineList.from_tumult_format(dump).reset_index()
                    if "strength" not in llist.columns:
                        llist["strength"] = 1.0
                    frame = llist[["lam0", "ion", "strength"]].copy()
                    frame["ion"] = frame.ion.astype(str).str.strip()
                except Exception:
                    logging.exception("Could not read the reference dump %s", dump)
            cache[dump] = frame
        return cache[dump]

    # --- the modal's run: reference, then corrections, as one chained job ---------------
    def run(self):
        obj = self._pending_atom
        if obj is None:
            return
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason or "Configure the target atmosphere.")
            pulse(self._modal, "reject")
            return
        ref_grid = self.grids.reference_grid()
        if ref_grid is None:
            self.notify("Flag a 1D LTE reference grid in Preferences → Atmosphere grids first.")
            pulse(self._modal, "reject")
            return
        # every option is captured now and carried through both stages: reopening the
        # modal while the reference runs must not alter the admitted calculation
        v = self.form1.values()
        opts = CorrectionOpts(profile=str(v["profile"]), ew_min=float(v["ew_min"]),
                              limit=float(v.get("limit") or DEFAULTS["limit"]),
                              line_grid=bool(v.get("line_grid", True)),
                              linelist=bool(v.get("linelist", False)))
        if opts.profile == "trapz" and not opts.line_grid:
            self.notify("The trapezoidal width needs the line grid for its windows; "
                        "switch it on or pick a fitted profile.")
            pulse(self._modal, "reject")
            return
        try:
            request = request.update(
                departures=DepartureRequest((obj.path,), line_grid=opts.line_grid))
        except ValueError as exc:
            self.notify(f"Cannot compute corrections: {exc}")
            pulse(self._modal, "reject")
            return
        self.draft.commit()
        pulse(self._modal, "close")
        self.notify(f"Computing the 1D LTE reference for {obj.symbol}…")
        # one durable record runs both stages; the stage-1 result arrives as its first
        # progress snapshot, every correction round as the next ones. The close patch above
        # only reaches the browser once this callback returns, so the actual submission - a
        # blocking round trip to the workflow daemon - is deferred one document tick (as
        # synthmodal.SynthModal._begin_submit/_end_submit does), rather than sitting in front
        # of it and leaving the modal looking inert for however long that round trip takes.
        kwargs = {"reference_atmosphere": atmosphere_grid(ref_grid),
                  "background_lines": opts.linelist, "profile": opts.profile,
                  "ew_min": opts.ew_min, "limit": opts.limit}
        name = self.ctx.unique_name(naming.corrections(obj.symbol, request,
                                                       self.setup.active_grid()))
        self.ctx.schedule(lambda: self.jobs.submit_operation(
            name, self.task_kind, request, lambda result: self._on_task_done(result, name),
            kwargs=kwargs, on_progress=lambda snapshot: self._on_task_progress(snapshot, name)))

    task_kind = "abundance.nlte_corrections"

    def _on_task_progress(self, snapshot, name):
        """A stage-1 result (the first snapshot of a run) is listed as ``name`` and its
        reference spectrum shown; a later snapshot updates the listed result in place."""
        if self.find(snapshot.uid) is None:
            self._on_reference_ready(snapshot, name)
        else:
            self._job_progress(self._adopt(snapshot, name))

    def _on_task_done(self, result, name):
        if self.find(result.uid) is None:
            self._on_reference_ready(result, name)  # its stage-1 snapshot never arrived
        self._on_corr_done(self._adopt(result, name))

    import_task_progress = lambda self, snapshot, handle: self._on_task_progress(
        snapshot, handle.label)
    import_task_result = lambda self, result, handle: self._on_task_done(result, handle.label)

    def _on_reference_ready(self, result, name):
        self._add_result(result, name)
        self.set_step(1)
        # draws the reference spectrum too, as this result's view is now selected
        self.select_result(result.uid)
        self.notify(f"{result.symbol}: {len(result.lines)} reference lines — computing "
                    "NLTE corrections…")

    def _on_corr_done(self, result):
        s = result.summary
        self.push_results()
        if self.selected() is result:
            self._render_selected()
            if self.jumped is not None and self.jumped < len(result.lines):
                self._render_line_decor(result.lines[self.jumped])
        if np.isfinite(s["a_mean"]):
            self.notify(f"{result.symbol}: ⟨ΔA⟩ = {s['a_mean']:+.3f} ± "
                        f"{self._summary_sigma(s):.3f} (n={s['n']}) — see the Scatter step.")

    # --- rendering ----------------------------------------------------------------------
    def _render_scatter(self):
        """Repaint both diagnostic panels (ΔA against E_low and against reduced EW) from
        the included, inverted lines; whiskers of zero length stand in for a missing
        formal error (legacy sessions)."""
        res = self.selected()
        inc = [(i, ln) for i, ln in enumerate(res.lines)
               if ln.include and np.isfinite(ln.a)] if res else []
        if not inc:
            self.panels.src.data = empty_points()
            self.panels.set_mean(None)
            return
        data = self.panels.pack_points(inc)
        self.panels.src.data = data
        self.panels.crop_panels((data,))
        s = res.summary
        self.panels.set_mean(s.get("a_mean", np.nan), self._summary_sigma(s))

    def table_html(self, res):
        columns = (
            ("λ₀ [Å]", "lam0"), ("ion", "ion"), ("E<sub>low</sub>", "e_low"),
            ("log gf", "loggf"), ("EW [mÅ]", "ew"), ("REW", "rew"), ("ΔA", "a"),
            ("flags", None))
        head = self._table_head(res, columns)
        body = []
        for i in self._line_order(res):
            ln = res.lines[i]
            chk = " checked" if ln.include else ""
            cls = "" if ln.include else " line-excluded"
            if self.jumped == i:
                cls += " row-selected"
            ew = format_value(ln.ew, ln.ew_err, digits=1)
            rew = format_value(reduced_ew(ln.ew, ln.lam0) if ln.ew > 0 else float("nan"))
            da = format_value(ln.a, ln.sigma_a if res.stage == "corr" else None, plus=True)
            fl = html.escape(", ".join(ln.flags))
            if window_edited(ln):
                fl += (f' <button class="text-btn" data-reset-window="{i}" title="Return to '
                       'the window this line was computed over">Reset window</button>')
            body.append(
                f'<tr class="corr-row{cls}" data-line="{i}">'
                f'<td><input type="checkbox" class="line-include-chk" data-line="{i}"{chk}></td>'
                f'<td>{ln.lam0:.3f}</td><td>{html.escape(ln.ion)}</td>'
                f'<td>{ln.e_low:.2f}</td><td>{ln.loggf:.2f}</td><td>{ew}</td>'
                f'<td>{rew}</td><td>{da}</td>'
                f'<td class="abund-flags">{fl}</td></tr>')
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
        return (f'<table class="linelist-table abund-table corr-lines" data-tool="corr">'
                f'<thead>{head}</thead>'
                f'<tbody>{"".join(body) or self._filtered_out_row(columns)}</tbody></table>{foot}')
