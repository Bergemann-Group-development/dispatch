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

The live objects, the focus and the collections' views are reached through the document's
:class:`.objects.WorkspaceObjects`; the filter tool and the display-settings form, whose
state a slot also records, are passed in beside it.
"""
from __future__ import annotations

from datetime import datetime
import html
import logging
import time

import numpy as np

from tumult.spectrum_io import Spectrum, load_file, load_run

from ..ui.latestwork import LatestWork
from ..ui.debounce import Trailing
from ..ui.htmlbridge import pulse
from . import session, store
from .objects import FOCUS_TOOLS
from .session import (AtomSnap, DepartureSnap, FilterSnap, FocusSnap, MaskSnap,
                      PreparedAtom, PreparedDepartures, PreparedSpectrum, PreparedWorkspace,
                      ResultsSnap, SpectrumSnap, ToolsSnap, ViewSnap, Workspace)

log = logging.getLogger("tumultstudio.jobs")

AUTOSAVE_DELAY_MS = 2000


def load_spectrum(obj):
    """A saved spectrum reference -> Spectrum (embedded synth data, a file, or a tumult run)."""
    if obj.kind == "synth":
        return Spectrum.from_dict(obj.data)
    return load_run(obj.path) if obj.kind == "run" else load_file(obj.path)


def load_atom(path):
    """A saved workspace's model atom (imported on demand: parsing an atom pulls in the
    whole atomic-data stack, which a document that opens none never needs)."""
    from tumult.atom import AtomClass
    return AtomClass(str(path))


def load_departures(folder):
    """A saved workspace's departure solve, reopened from its run folder. The populations
    stay on disk; only the record is read here."""
    from tumult.departures import Departures
    return Departures.load(str(folder))


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

    ``workspace`` is the document's :class:`.objects.WorkspaceObjects`; ``filter`` the
    filter tool and ``settings`` the display-settings form, whose state a slot records too;
    ``refresh`` the document's full repaint after a restore. ``deserialise_results`` and
    ``relink_results`` build and relink a saved result section (:mod:`..analysis.results`);
    with :func:`load_spectrum` and :func:`load_atom` they are plain functions used off the
    document thread, so a slot is interpreted without constructing any view.
    """

    def __init__(self, ctx, workspace, delivery, workflows, *, filter, settings, refresh,
                 deserialise_results, relink_results,
                 load_spectrum=load_spectrum, load_atom=load_atom,
                 load_departures=load_departures):
        self.ctx = ctx
        self.workspace = workspace
        self.delivery = delivery
        self.workflows = workflows
        self.filter = filter
        self.settings = settings
        self.refresh = refresh
        self._load_spectrum = load_spectrum
        self._load_atom = load_atom
        self._load_departures = load_departures
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
        ws, filtertool = self.workspace, self.filter
        canvas = ws.canvas
        spectra, atomtool, deptool = ws.spectra, ws.atoms, ws.departures
        specs = []
        for o in spectra:
            kind, path = o.source if o.source else (None, None)
            if kind == "synth":
                data, path = o.raw.to_dict(), ""
            elif kind in ("file", "run") and path:
                data = None
            else:
                continue          # nothing to restore it from
            knots = ({"x": list(map(float, o.norm_knots[0])),
                      "y": list(map(float, o.norm_knots[1]))}
                     if o.norm_knots is not None else None)
            # synth params are tiny and have no source file, so persist them; FITS metadata
            # is large and re-read from disk on restore, so isn't stored.
            specs.append(SpectrumSnap(
                uid=o.uid, kind=kind, name=o.name, path=path, data=data,
                meta=o.meta if kind == "synth" else None,
                color=o.color, dash=o.dash, width=o.width, alpha=o.alpha,
                visible=o.visible, conv=dict(o.conv), norm_knots=knots,
                already_normalised=bool(o.already_normalised)))
        saved_uids = {sp.uid for sp in specs}
        spectrum_selected = (spectra.selected_uid if spectra.selected_uid in saved_uids
                             else None)
        configured = filtertool._configured
        # only persist a low-pass cutoff when a low-pass filter actually exists
        has_lowpass = (filtertool.method == "lowpass"
                       and any(o.filtered is not None for o in spectra))
        departures_snap = deptool.snapshot()
        focus_kind, focus_uid = ws.focused()
        if focus_kind == "spectrum":
            focus_uid = spectrum_selected
        wizards = ws.results
        return Workspace(
            spectra=specs, spectrum_selected=spectrum_selected,
            masks=MaskSnap(**ws.masks.snapshot()),
            atoms=[AtomSnap(uid=a.uid, path=a.path) for a in atomtool.atoms],
            atom_selected=atomtool.selected_uid,
            departures=[DepartureSnap(**d) for d in departures_snap["departures"]],
            departure_selected=departures_snap["selected"],
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
                               load_departures=self._load_departures,
                               deserialise=self._deserialise_results,
                               relink=self._relink_results)

    def capture_live(self):
        """The live workspace as a prepared one (no file reads): what a failed commit
        re-applies."""
        saved = Workspace.from_dict(self.snapshot())
        ws = self.workspace
        by_uid = {o.uid: o for o in ws.spectra}
        return PreparedWorkspace(
            spectra=[PreparedSpectrum(sp, by_uid[sp.uid].raw, by_uid[sp.uid].meta)
                     for sp in saved.spectra],
            spectrum_selected=saved.spectrum_selected,
            mask_sets=list(ws.masks.sets), mask_current=ws.masks.current_uid,
            mask_active=ws.masks.active,
            atoms=[PreparedAtom(AtomSnap(uid=a.uid, path=a.path), a.atom)
                   for a in ws.atoms.atoms],
            atom_selected=ws.atoms.selected_uid,
            departures=[PreparedDepartures(snap, ws.departures.find(snap.uid).departures)
                        for snap in saved.departures],
            departure_selected=saved.departure_selected,
            results={name: (list(w.results), w.results.selected)
                     for name, w in ws.results.items()},
            view=saved.view, filter=saved.filter, focus=saved.focus, tools=saved.tools)

    def _place_spectrum(self, ps):
        canvas = self.workspace.canvas
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
        ctx, ws, filtertool = self.ctx, self.workspace, self.filter
        canvas, spectra, wizards = ws.canvas, ws.spectra, ws.results
        ctx.tools.deactivate_active()       # whatever view was open, the canvas returns
        for w in wizards.values():
            w.blur()
        ws.focus_kind = "spectrum"
        filtertool.clear()
        with canvas.batch_update():
            for o in list(spectra):
                canvas.remove_object(o.uid)
            for ps in p.spectra:
                self._place_spectrum(ps)
        spectra.selected_uid = p.spectrum_selected     # None when none was, or it is gone
        ws.masks.restore(p.mask_sets, p.mask_current, p.mask_active)
        ws.atoms.restore(p.atoms, p.atom_selected, p.tools.atom_tab)
        ws.departures.restore(p.departures, p.departure_selected)
        for name, w in wizards.items():
            results, selected = p.results[name]
            w.restore(results, selected)
            if name in p.tools.steps:
                w.set_step(p.tools.steps[name])
        ws.atoms.render_objects()
        v = p.view
        if v.xunit:
            canvas.set_xunit(v.xunit)
        if v.yunit:
            canvas.set_yunit(v.yunit)
        self.settings.set_values(xunit=canvas.xlabel, yunit=canvas.ylabel)  # sync the form
        f = p.filter
        if spectra.objects and f.method == "lowpass" and f.res is not None:
            filtertool.apply_saved_lowpass(f.res)
            if not f.active:
                filtertool.set_visible(False)
        elif spectra.objects and f.fmax is not None and f.floor is not None:
            filtertool.apply_saved(f.fmax, f.floor)
            if not f.active:
                filtertool.set_visible(False)
        # the focus, through the same routing a click takes, so the panels it opens are
        # the ones the user had; a cleared focus leaves the plain canvas
        ws.focus_object(p.focus.kind, p.focus.uid)
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
        if active in FOCUS_TOOLS and ctx.tools.active_name() != active:
            ctx.tools.activate(active)
        # the explicit ranges last, after everything above that moves the view
        _restore_range(canvas.nav.x_range, v.nav_x)
        _restore_range(canvas.main.x_range, v.main_x)
        _restore_range(canvas.main.y_range, v.main_y)
        _restore_range(canvas.nav.y_range, v.nav_y)
        self.refresh()

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
        self.refresh()
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
