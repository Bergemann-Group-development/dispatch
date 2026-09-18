"""Saving and restoring one document's workspace.

:mod:`.session` defines what a saved workspace *is* — the JSON schema, its validation and
the staged preparation that turns a slot into loaded objects — and :mod:`.store` owns the
config file it lives in. This module is the document's side of that: it takes the
snapshot of the live collections, drives the validate → prepare → commit flow with a
rollback to the workspace that was there, reports what a slot could not restore, and
keeps the autosave slot current so a closed tab resumes.

A commit is a workspace boundary. The generation advances, so anything still running for
the old workspace settles its row and is dropped rather than landing in the new one; this
document's own non-detached jobs are cancelled; every collection is replaced, empty ones
included; and the slot's durable-workflow references replace this document's, which is
where reattachment happens (:mod:`..execution.workflows`). If applying the prepared
workspace fails, the live one captured beforehand is put back and the failure is
reported; if that fails too, the document says so rather than pretending.

The live collections are reached through :class:`WorkspaceParts`, the seam where a single
owner of live objects and focus belongs. Until that owner exists, this is the explicit
list of what a workspace is made of, rather than a closure over the whole document.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import html
import logging
import time
from typing import Callable

import numpy as np

from ..execution.latestwork import LatestWork
from ..ui.debounce import Trailing
from ..ui.htmlbridge import pulse
from . import session, store
from .session import (AtomSnap, FilterSnap, FOCUS_TOOL, FocusSnap, MaskSnap,
                      PreparedAtom, PreparedSpectrum, PreparedWorkspace, ResultsSnap,
                      SpectrumSnap, ToolsSnap, ViewSnap, Workspace)

log = logging.getLogger("tumultstudio.jobs")

AUTOSAVE_DELAY_MS = 2000


@dataclass(frozen=True)
class WorkspaceParts:
    """The live collections of one document that a saved workspace covers.

    The canvas, the tool manager and the mask store are reached through the context;
    these are the rest. ``focused`` reports the single left-sidebar focus as
    ``(kind, uid)``, ``focus_spectrum`` routes a focus change the way a click does, and
    ``refresh`` is the document's full repaint after a restore.
    """
    atoms: object                        # AtomTool
    mask_view: object                    # MaskTool (the line-mask view)
    results: dict                        # tool name -> the result view owning that list
    filter: object                       # FilterTool
    settings: object                     # the display-settings Form
    focused: Callable[[], tuple]
    focus_spectrum: Callable[[str], None]
    refresh: Callable[[], None]


def _range(r):
    """``[start, end]`` of a Bokeh range when it holds an explicit finite interval."""
    lo, hi = r.start, r.end
    if not all(isinstance(v, (int, float)) and np.isfinite(v) for v in (lo, hi)):
        return None
    return [float(lo), float(hi)] if lo < hi else None


def _restore_range(rng, vals):
    if not vals:
        return
    lo, hi = vals
    if rng.bounds:
        lo, hi = max(lo, min(rng.bounds)), min(hi, max(rng.bounds))
    if lo < hi:
        rng.start, rng.end = lo, hi


def sessions_html(slots):
    """Render the saved-session slots as a clickable list (game-style save menu)."""
    if not slots:
        return '<div class="file-empty">No saved sessions yet — make one below.</div>'
    rows = []
    for s in slots:
        try:
            when = datetime.fromtimestamp(s.get("saved_at", 0)).strftime("%Y-%m-%d %H:%M")
        except (OSError, ValueError, OverflowError):
            when = ""
        sub = session.summary(s.get("snapshot")) + (f' · {when}' if when else "")
        auto = s.get("autosave", False)     # the running document's own slot, listed first
        rows.append(
            f'<div class="session-item{" autosave" if auto else ""}" data-sid="{s["id"]}" '
            f'tabindex="0" role="button" title="{"The last session, kept current while a document runs" if auto else "Select saved workspace"}">'
            f'<i class="ti ti-{"history" if auto else "device-floppy"}"></i>'
            f'<div class="session-meta"><span class="session-name">{html.escape(s["name"])}</span>'
            f'<span class="session-sub">{html.escape(sub)}</span></div>'
            f'<button class="session-del" data-sid="{s["id"]}" '
            f'title="Delete saved workspace"><i class="ti ti-trash"></i></button>'
            f'</div>')
    return '<div class="session-list">' + "".join(rows) + "</div>"


class WorkspacePersistence:
    """Snapshot, restore and autosave of one document's workspace.

    ``load_spectrum`` reads a saved spectrum reference (embedded synthesis data, a file or
    a tumult run) and ``load_atom`` a model atom; ``deserialise_results`` and
    ``relink_results`` build and relink a saved result section (:mod:`..analysis.results`).
    All four are plain functions used off the document thread, so a slot is interpreted
    without constructing any view.
    """

    def __init__(self, ctx, parts, delivery, workflows, *, load_spectrum, load_atom,
                 deserialise_results, relink_results):
        self.ctx = ctx
        self.parts = parts
        self.delivery = delivery
        self.workflows = workflows
        self._load_spectrum = load_spectrum
        self._load_atom = load_atom
        self._deserialise_results = deserialise_results
        self._relink_results = relink_results

        # a dedicated space-delimited channel: names may contain spaces/colons, so the
        # ":"-delimited object-actions bridge can't carry them safely (cf. cache-actions).
        self._actions = ctx.bridges.action("session-actions")
        # broadcast to every session list on the page: the Sessions modal *and* the inline
        # empty-state loader shown over the canvas when nothing is loaded
        self._sessions = ctx.bridges.html_all(".js-session-results", name="session-render")
        self._sessions.value = sessions_html(store.list_sessions())  # seed initial list
        self._load_prompt = ctx.bridges.dialog("workspace-load-modal")
        self._sessions_dialog = ctx.bridges.dialog("session-modal")
        self._actions.on_change("value", self._on_action)

        self._saved = {"snapshot": self.snapshot(), "pending": None}

        # --- autosave: the live workspace kept in one slot, so a closed tab resumes -----
        # Every refresh (a load, a result, a view change) schedules a write once the burst
        # settles; the session-destroyed hook writes at once. A blank workspace is never
        # written, so a fresh tab does not replace the last real session, and a snapshot
        # equal to the last one written is skipped. Failures are logged: this must never
        # interrupt the operation that triggered it.
        self._autosave_delay = Trailing(AUTOSAVE_DELAY_MS)
        self._autosaved = {"snapshot": None, "requested": None}
        self._autosave_work = LatestWork(
            lambda *args: ctx.submit_io(*args), "autosave workspace",
            self._autosave_failed,
            immediate=getattr(ctx.doc, "session_context", None) is None)

        ctx.restore_workspace = self.restore   # local harnesses load slots through it
        ctx.workspace_changed = self.autosave

    @property
    def _headless(self):
        """True in a local harness: there is no client session, so the work a Bokeh
        server would hand to a pool runs inline and the caller sees the outcome."""
        return getattr(self.ctx.doc, "session_context", None) is None

    # --- the snapshot ----------------------------------------------------------------
    def snapshot(self):
        """The live workspace as a :class:`session.Workspace` dict."""
        canvas, parts = self.ctx.canvas, self.parts
        filtertool, atomtool = parts.filter, parts.atoms
        specs = []
        for o in canvas.objects:
            kind, path = o.source if o.source else (None, None)
            if kind == "synth":
                data, path = o.raw.to_dict(), ""
            elif kind in ("file", "run") and path:
                data = None
            else:
                continue          # nothing to restore it from
            g = o.line.glyph
            knots = ({"x": list(map(float, o.norm_knots[0])),
                      "y": list(map(float, o.norm_knots[1]))}
                     if o.norm_knots is not None else None)
            # synth params are tiny and have no source file, so persist them; FITS metadata
            # is large and re-read from disk on restore, so isn't stored.
            specs.append(SpectrumSnap(
                uid=o.uid, kind=kind, name=o.name, path=path, data=data,
                meta=o.meta if kind == "synth" else None,
                color=g.line_color, dash=o.dash, width=g.line_width, alpha=g.line_alpha,
                visible=o.visible, conv=dict(o.conv), norm_knots=knots,
                already_normalised=bool(o.already_normalised)))
        saved_uids = {sp.uid for sp in specs}
        spectrum_selected = canvas.selected_uid if canvas.selected_uid in saved_uids else None
        configured = filtertool._configured
        # only persist a low-pass cutoff when a low-pass filter actually exists
        has_lowpass = (filtertool.method == "lowpass"
                       and any(o.filtered is not None for o in canvas.objects))
        focus_kind, focus_uid = parts.focused()
        if focus_kind == "spectrum":
            focus_uid = spectrum_selected
        wizards = parts.results
        return Workspace(
            spectra=specs, spectrum_selected=spectrum_selected,
            masks=MaskSnap(**self.ctx.masks.snapshot()),
            atoms=[AtomSnap(uid=a.uid, path=a.path) for a in atomtool.atoms],
            atom_selected=atomtool.selected_uid,
            results={name: ResultsSnap(**w.snapshot()) for name, w in wizards.items()},
            view=ViewSnap(xunit=canvas.xlabel, yunit=canvas.ylabel,
                          main_x=_range(canvas.main.x_range), main_y=_range(canvas.main.y_range),
                          nav_x=_range(canvas.nav.x_range), nav_y=_range(canvas.nav.y_range)),
            filter=FilterSnap(
                active=canvas.show_filtered, method=filtertool.method,
                res=filtertool.form.value("res") if has_lowpass else None,
                fmax=filtertool.fmax if configured else None,
                floor=filtertool.floor if configured else None),
            focus=FocusSnap(kind=focus_kind, uid=focus_uid),
            tools=ToolsSnap(active=self.ctx.tools.active_name(), atom_tab=atomtool.tab,
                            steps={name: w.step for name, w in wizards.items()},
                            jumped={name: w.jumped for name, w in wizards.items()}),
        ).to_dict()

    # --- validate, prepare, commit ----------------------------------------------------
    def prepare(self, snap):
        """Validate a slot and build its objects off the live workspace
        (:func:`session.prepare`); :class:`session.SessionError` when it is unusable."""
        ws = Workspace.from_dict(snap)
        return session.prepare(ws, load_spectrum=self._load_spectrum,
                               load_atom=self._load_atom,
                               deserialise=self._deserialise_results,
                               relink=self._relink_results)

    def capture_live(self):
        """The live workspace as a prepared one (no file reads): what a failed commit
        re-applies."""
        ws = Workspace.from_dict(self.snapshot())
        by_uid = {o.uid: o for o in self.ctx.canvas.objects}
        wizards = self.parts.results
        return PreparedWorkspace(
            spectra=[PreparedSpectrum(sp, by_uid[sp.uid].raw, by_uid[sp.uid].meta)
                     for sp in ws.spectra],
            spectrum_selected=ws.spectrum_selected,
            mask_sets=list(self.ctx.masks.sets), mask_current=self.ctx.masks.current_uid,
            mask_active=self.ctx.masks.active,
            atoms=[PreparedAtom(AtomSnap(uid=a.uid, path=a.path), a.atom)
                   for a in self.parts.atoms.atoms],
            atom_selected=self.parts.atoms.selected_uid,
            results={name: (list(w.results), w.results.selected)
                     for name, w in wizards.items()},
            view=ws.view, filter=ws.filter, focus=ws.focus, tools=ws.tools)

    def _place_spectrum(self, ps):
        canvas = self.ctx.canvas
        snap, spec = ps.snap, ps.spectrum
        o = canvas.add_object(spec, name=snap.name, color=snap.color, dash=snap.dash,
                              source=(snap.kind, snap.path or None), meta=ps.meta,
                              uid=snap.uid)
        if snap.width is not None or snap.alpha is not None:
            canvas.set_style(o, width=snap.width, alpha=snap.alpha)
        knots = snap.norm_knots
        o.already_normalised = snap.already_normalised
        if knots and knots["x"]:
            o.norm_knots = (np.asarray(knots["x"], float), np.asarray(knots["y"], float))
        if o.already_normalised or (knots and knots["x"]):
            canvas.normalise_object(o)
        o.conv = dict(snap.conv)
        if o.conv or o.already_normalised or (knots and knots["x"]):
            canvas._refresh_object(o)  # apply conv/normalisation to main + nav
        if not snap.visible:
            canvas.toggle_object(o.uid)

    def _apply(self, p):
        """Replace every collection with the prepared workspace ``p`` and restore its
        view, selections and focus. Every collection is replaced, empty ones included;
        spectra, masks and atoms land before the results that refer to them."""
        ctx, parts = self.ctx, self.parts
        canvas, wizards = ctx.canvas, parts.results
        ctx.tools.deactivate_active()       # whatever view was open, the canvas returns
        for w in wizards.values():
            w.blur()
        ctx.focus_kind = "spectrum"
        parts.filter.clear()
        with canvas.batch_update():
            for o in list(canvas.objects):
                canvas.remove_object(o.uid)
            for ps in p.spectra:
                self._place_spectrum(ps)
        canvas.selected_uid = p.spectrum_selected     # None when none was, or it is gone
        ctx.masks.restore(p.mask_sets, p.mask_current, p.mask_active)
        parts.atoms.restore(p.atoms, p.atom_selected, p.tools.atom_tab)
        for name, w in wizards.items():
            results, selected = p.results[name]
            w.restore(results, selected)
            if name in p.tools.steps:
                w.set_step(p.tools.steps[name])
        parts.atoms.render_objects()
        v = p.view
        if v.xunit:
            canvas.set_xunit(v.xunit)
        if v.yunit:
            canvas.set_yunit(v.yunit)
        parts.settings.set_values(xunit=canvas.xlabel, yunit=canvas.ylabel)  # sync the form
        f = p.filter
        if canvas.objects and f.method == "lowpass" and f.res is not None:
            parts.filter.apply_saved_lowpass(f.res)
            if not f.active:
                parts.filter.set_visible(False)
        elif canvas.objects and f.fmax is not None and f.floor is not None:
            parts.filter.apply_saved(f.fmax, f.floor)
            if not f.active:
                parts.filter.set_visible(False)
        # the focus, through the same routing a click takes, so the panels it opens are
        # the ones the user had; a cleared focus leaves the plain canvas
        focus = p.focus
        if focus.kind == "spectrum":
            if focus.uid is not None:
                parts.focus_spectrum(focus.uid)
        elif focus.kind == "linemask":
            parts.mask_view.select_set(focus.uid)
        elif focus.kind == "atom":
            parts.atoms.select(focus.uid)
        elif focus.uid is not None:
            wizards[FOCUS_TOOL[focus.kind]].select_result(focus.uid)
        for name, w in wizards.items():
            jumped = p.tools.jumped.get(name)
            res = w.selected()
            if jumped is None or res is None or jumped >= len(res.lines):
                continue
            if w.active and w.is_focused():
                w.jump_line(str(jumped))       # the accent, hatch and CoG of that line
            else:
                w.jumped = jumped              # painted when its view is next opened
        # a panel open beside the canvas (a wizard's steps, the mask or atom view) is
        # reopened; the range, filter and normalise gestures are transient and stay closed
        active = p.tools.active
        if active in ("abundance", "abundcorr", "ionbal", "linemask", "atom") \
                and ctx.tools.active_name() != active:
            ctx.tools.activate(active)
        # the explicit ranges last, after everything above that moves the view
        _restore_range(canvas.nav.x_range, v.nav_x)
        _restore_range(canvas.main.x_range, v.main_x)
        _restore_range(canvas.main.y_range, v.main_y)
        _restore_range(canvas.nav.y_range, v.nav_y)
        parts.refresh()

    def commit(self, prepared, name="save", workflows=(), live=None):
        """Commit an already-prepared slot on the document thread with rollback."""
        ctx = self.ctx
        if live is None:
            live = self.capture_live()
        ctx.advance_generation()
        cancelled = [job.name for job in ctx.jobs.active if not job.detached]
        self.delivery.cancel_document_work()
        try:
            self._apply(prepared)
        except Exception as exc:
            logging.exception("Applying the saved workspace %r failed; restoring the "
                              "previous workspace", name)
            try:
                self._apply(live)
            except Exception:
                logging.exception("The previous workspace could not be restored either")
                ctx.notify(f"Loading “{name}” failed: {exc}; the previous workspace could "
                           "not be restored", level="danger", revert_ms=10000)
                return None
            ctx.notify(f"Loading “{name}” failed: {exc}; the previous workspace was restored",
                       level="danger", revert_ms=10000)
            return None
        outcomes, gone = self.workflows.reattach(workflows)
        report = session.LoadReport(name, prepared.counts(), prepared.skipped + gone,
                                    prepared.unresolved, cancelled, outcomes)
        self.parts.refresh()
        for line in report.details():
            log.warning("%s", line)
        ctx.notify(report.message(), level="warn" if report.partial else "success",
                   revert_ms=12000 if report.partial else 4000)
        return report

    def restore(self, snap, name="save", workflows=()):
        """Load a slot synchronously (used by local harnesses).

        The interactive session action prepares on a worker; both paths use the same
        validated commit and rollback. The commit marks a workspace boundary (new
        generation, old result subscriptions forgotten, this document's non-detached jobs
        cancelled) and replaces every collection. If applying fails, the previous
        workspace is put back. The slot's ``workflows`` references then replace this
        document's (see :meth:`..execution.workflows.WorkflowService.reattach`). Returns
        the :class:`session.LoadReport`, or ``None`` when nothing was loaded."""
        try:
            prepared = self.prepare(snap)
            live = self.capture_live()
        except Exception as exc:   # SessionError names the field; anything else is logged
            logging.exception("The saved workspace %r cannot be loaded", name)
            self.ctx.notify(f"Could not load “{name}”: {exc}", level="danger", revert_ms=10000)
            return None
        return self.commit(prepared, name, workflows, live)

    def restore_async(self, snap, name="save", workflows=(), on_done=None):
        """Prepare a saved workspace on a worker, then commit the current generation."""
        try:
            live = self.capture_live()
        except Exception as exc:
            logging.exception("The current workspace could not be captured")
            self.ctx.notify(f"Could not start loading “{name}”: {exc}", level="danger")
            return None

        def prepared(value):
            report = self.commit(value, name, workflows, live)
            if report is not None and on_done is not None:
                on_done(report)

        return self.delivery.load_async(f"load {name}", lambda token: self.prepare(snap),
                                        prepared)

    # --- autosave ---------------------------------------------------------------------
    def _autosave_failed(self, exc):
        self._autosaved["requested"] = None
        logging.error("Autosave failed: %s", exc)

    def _write_autosave(self, synchronous=False):
        try:
            snap = self.snapshot()
            if (session.is_blank(snap) or snap == self._autosaved["snapshot"]
                    or snap == self._autosaved["requested"]):
                return
            workflows = self.ctx.workflow_refs.to_list()
            if synchronous:
                store.save_autosave(snap, workflows=workflows)
                self._autosaved["snapshot"] = snap
                return
            self._autosaved["requested"] = snap

            def work(token):
                store.save_autosave(snap, workflows=workflows)
                return snap

            def done(saved):
                self._autosaved["snapshot"] = saved
                if self._autosaved["requested"] == saved:
                    self._autosaved["requested"] = None

            self._autosave_work.request(work, done)
        except Exception:
            logging.exception("Autosave failed")

    def autosave(self, now=False):
        """Schedule the autosave write once the current burst of changes settles, or
        write it at once (the document is closing)."""
        if now:
            self._autosave_delay.cancel("autosave")
            self._write_autosave(synchronous=True)
        else:
            self._autosave_delay("autosave", self._write_autosave)

    # --- the saved slots ---------------------------------------------------------------
    def slots_html(self):
        """The saved slots as the page first paints them."""
        return sessions_html(store.list_sessions())

    def render_sessions(self):
        """Re-read the slots into every session list; a config the store rejects leaves
        the list as it was and reports why."""
        try:
            self._sessions.value = sessions_html(store.list_sessions())
        except store.ConfigError as exc:
            logging.exception("Could not list the saved sessions")
            self.ctx.notify(f"Could not read saved sessions: {exc}", level="danger")

    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]  # strip the trailing nonce
        # every store call may raise ConfigError (unreadable/invalid config, failed
        # replacement, lock timeout); the message names the path and the cause, and the
        # live workspace stays as it is so the user can retry after fixing the cause.
        if kind == "save":
            self._save_slot(arg)
        elif kind in ("load", "load-discard", "load-save"):
            self._load_slot(kind, arg)
        elif kind == "delete":
            self._delete_slot(arg)
        elif kind == "refresh":
            self.render_sessions()  # re-read slots (another tab may have changed them)

    def _save_slot(self, arg):
        name = arg.strip() or time.strftime("Save %Y-%m-%d %H:%M")
        snap = self.snapshot()
        workflows = self.ctx.workflow_refs.to_list()

        def saved(slot):
            self._saved["snapshot"] = snap
            self.render_sessions()
            self.ctx.notify(f"Saved “{name}”", level="success")

        if self._headless:
            saved(store.save_session(name, snap, workflows=workflows))
        else:
            self.ctx.submit_io(f"save {name}",
                               lambda token: store.save_session(name, snap,
                                                                workflows=workflows),
                               saved)

    def _load_slot(self, kind, arg):
        if kind == "load" and self.snapshot() != self._saved["snapshot"]:
            self._saved["pending"] = arg
            pulse(self._load_prompt, "open")
            return
        if kind != "load":
            arg = self._saved["pending"]
            if arg is None:
                return
            if kind == "load-save":
                try:
                    store.save_session(time.strftime("Workspace %Y-%m-%d %H:%M"),
                                       self.snapshot(),
                                       workflows=self.ctx.workflow_refs.to_list())
                except store.ConfigError as exc:
                    self.ctx.notify(f"Could not save workspace: {exc}", level="danger")
                    return
            pulse(self._load_prompt, "close")
        try:
            slot = store.get_slot(int(arg))
        except store.ConfigError as exc:
            logging.exception("Could not read the saved session %s", arg)
            self.ctx.notify(f"Could not read the save: {exc}", level="danger")
            return
        if slot is None:
            self.ctx.notify("That save is no longer available.", level="danger")
            self.render_sessions()
            return

        def loaded(report):
            self._saved["snapshot"] = self.snapshot()
            self._saved["pending"] = None
            pulse(self._sessions_dialog, "close")

        if self._headless:
            report = self.restore(slot["snapshot"], slot["name"], slot["workflows"])
            if report is not None:
                loaded(report)
        else:
            self.restore_async(slot["snapshot"], slot["name"], slot["workflows"], loaded)

    def _delete_slot(self, arg):
        try:
            found = store.delete_session(int(arg))
        except store.ConfigError as exc:
            logging.exception("Could not delete the saved session %s", arg)
            self.ctx.notify(f"Delete failed: {exc}", level="danger")
            return
        self.render_sessions()
        if not found:
            self.ctx.notify("That save was already gone.", level="warn")
