"""TumultStudio2 document construction, callable by Bokeh and local harnesses.

Run with:   bokeh serve --show TumultStudio2
"""
import html
import itertools
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

_nonce_counter = itertools.count()


def _nonce():
    return next(_nonce_counter)

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))  # repo root → `import tumult` (during startup)

import astropy.units as u  # noqa: E402
import numpy as np  # noqa: E402
from bokeh.io import curdoc  # noqa: E402
from bokeh.layouts import column  # noqa: E402

from tumult.spectrum_io import Spectrum, load_spectrum, load_file, is_run, load_run  # noqa: E402
from canvas import SpectrumCanvas, DASH_CSS, LINE  # noqa: E402
from context import AppContext  # noqa: E402
from controls import Field, Form, FORM_SCRIPT, NUM_INPUT_SCRIPT  # noqa: E402
from inspector import Inspector, swatch_html  # noqa: E402
from filtertool import FilterTool  # noqa: E402
from normalisetool import NormaliseTool  # noqa: E402
from maskstore import MaskStore  # noqa: E402
from masktool import MaskTool  # noqa: E402
from linenav import LineNavigator  # noqa: E402
from lineannot import LineAnnotations  # noqa: E402
from tumult import refstars  # noqa: E402
from tumult import linestrength  # noqa: E402
from tumult.linestrength import LineStrengths  # noqa: E402
from rangetool import RangeTool  # noqa: E402
from refstarspage import RefStarsPage  # noqa: E402
from synthmodal import SynthModal  # noqa: E402
from synthsetup import ATMOS_SCRIPT, SynthSetup  # noqa: E402
from abundancetool import AbundanceTool  # noqa: E402
from abundwizard import ANALYSIS_SCRIPT, PROFILES  # noqa: E402  (shared EW-profile choices)
from atomtool import AtomTool, AtomObject  # noqa: E402
from abundcorrtool import AbundCorrTool  # noqa: E402
from ionbalance import IonBalanceTool  # noqa: E402
from tumult.synthesis import workflows as synth  # noqa: E402  (CancelToken/CancelledError for cancellable loads)
from tumult.synthesis import submit as submit_synthesis  # noqa: E402
from tumult.synthesis import submit_task as submit_task_record  # noqa: E402
from tumult.synthesis import submission  # noqa: E402  (open workflows on reconnect)
from tumult.synthesis import daemon as workflow_daemon  # noqa: E402
from synthprefs import PREF_LABELS, load_prefs, _resolve_path  # noqa: E402  (the configured workflow root)
from filepicker import FilePicker  # noqa: E402
from toolmanager import ToolManager  # noqa: E402
import session  # noqa: E402
from session import (Workspace, SpectrumSnap, MaskSnap, AtomSnap, ResultsSnap, ViewSnap,  # noqa: E402
                     FilterSnap, FocusSnap, ToolsSnap, PreparedWorkspace, PreparedSpectrum,
                     PreparedAtom, FOCUS_TOOL, RESULT_TOOLS)
from htmlbridge import BridgeSet, pulse, toggle  # noqa: E402
from jobs import (ExternalRun, Job, Recoverable, QUEUED, RUNNING, CANCELLING, FINISHED, FAILED,  # noqa: E402
                  CANCELLED)
import processes  # noqa: E402
from processes import processes_html  # noqa: E402
from debounce import Trailing  # noqa: E402
from units import X_UNITS, Y_UNITS, to_range  # noqa: E402
from tumult import linestrength
from tumult.synthesis.runtime import (IO as _EXEC, WORKFLOWS as _SYNTH_EXEC, HOP as _HOP_EXEC,
                                      get_runtime as _get_runtime)
from profiling import profiled  # noqa: E402
import store  # noqa: E402

# how long the view has to hold still before the visible band's line strengths are computed
STRENGTH_DELAY_MS = 400

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(name)s: %(message)s")
# Job failures are logged here from the worker thread, traceback included, whether or not
# the document that admitted them still exists; the user-facing report goes through the
# owning document's Processes row and footer alone. No root-logger handler is installed
# per document, so one session's errors never reach another session's footer.
log = logging.getLogger("tumultstudio.jobs")


def picker_start_dir():
    """Where the file picker opens (the canvas itself starts empty): the folder of
    $TUMULT_SPECTRUM if set, else the most recently opened file's folder, else home."""
    env = os.environ.get("TUMULT_SPECTRUM")
    for candidate in ([env] if env else []) + store.recent_files():
        folder = Path(candidate).parent
        if folder.is_dir():
            return folder
    return Path.home()


def workflow_roots():
    """Where durable workflow records are looked for: the synthesis output root the
    preferences configure (``<dispatch_cwd>/<synth_output_dir>``), which is where
    ``tumult.synthesis.submit`` writes ``_workflows/``. Discovery never leaves it."""
    prefs = load_prefs()
    return [synth.execution_root(prefs) / (prefs.get("synth_output_dir") or "synth")]


def _load_obj(obj):
    """Load a session object reference -> Spectrum (embedded synth data, file, or tumult run)."""
    if obj.kind == "synth":
        return Spectrum.from_dict(obj.data)
    return load_run(obj.path) if obj.kind == "run" else load_file(obj.path)


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


def build_settings(canvas, on_units=None):
    def on_change(key, value):
        if key == "xunit":
            canvas.set_xunit(value)
        if on_units:
            on_units()

    # the y axis is always the normalised (dimensionless) flux — every spectrum is
    # continuum-normalised on load — so only the x unit is selectable
    return Form("settings", [
        Field("xunit", "X axis unit", kind="select",
              value=canvas.xlabel, options=list(X_UNITS)),
    ], on_change=on_change)


def objects_html(canvas, active_kind="spectrum"):
    rows = []
    for obj in canvas.objects:
        style = DASH_CSS.get(obj.dash, "dashed")
        color = obj.line.glyph.line_color if obj.line else LINE
        name = html.escape(obj.name)
        eye = "eye" if obj.visible else "eye-off"
        dim = " hidden-object" if not obj.visible else ""
        # only one item across the Objects + Abundances lists is highlighted at a time
        sel = " selected" if (obj.uid == canvas.selected_uid
                              and active_kind == "spectrum") else ""
        rows.append(
            f'<div class="object-item{dim}{sel}" data-uid="{obj.uid}">'
            f'<span class="dash-swatch" style="border-top-style:{style};border-top-color:{color}"></span>'
            f'<span class="object-name" title="{name}">{name}</span>'
            f'<button class="object-btn" data-action="toggle" data-uid="{obj.uid}" '
            f'title="Hide / show"><i class="ti ti-{eye}"></i></button>'
            f'<button class="object-btn" data-action="remove" data-uid="{obj.uid}" '
            f'title="Remove from workspace (the file stays on disk)"><i class="ti ti-x"></i></button>'
            f'</div>'
        )
    return "".join(rows)


# shown in the Objects panel while it holds nothing (and nothing is loading): how to
# get a first spectrum onto the canvas
OBJECTS_EMPTY_HINT = (
    '<div class="object-empty">No spectra yet. Open one with '
    '<b>Open → Spectrum</b> <i class="ti ti-folder-open"></i> above, or '
    '<b>New → Synthesise spectrum</b> to compute one.</div>'
)


def objects_panel_html(canvas, active_kind="spectrum"):
    """The Objects panel body: the loaded objects, or the empty-state hint when there is
    nothing to show. Work still running lives in the Processes panel (see
    :func:`processes_html`), so a burst of background tasks cannot bury the object list."""
    return objects_html(canvas, active_kind) or OBJECTS_EMPTY_HINT


def _field(label, value):
    v = html.escape(str(value))
    return (f'<div class="field-group"><label>{label}</label>'
            f'<div class="field-static" title="{v}">{v}</div></div>')


def _meta_html(meta):
    """Render an object's Info-tab provenance: a synth-parameter list or a scrollable
    table of FITS header cards (or nothing when there's no metadata)."""
    if not meta:
        return ""
    if meta.get("type") == "synth":
        rows = "".join(_field(label, value) for label, value in meta.get("params", []))
        return '<div class="section-header">Synthesis</div>' + rows
    if meta.get("type") == "fits":
        cards = meta.get("cards") or []
        if not cards:
            return ""
        body = "".join(
            f'<tr><td class="meta-key">{html.escape(k)}</td>'
            f'<td class="meta-val" title="{html.escape(v)}">{html.escape(v)}</td></tr>'
            for k, v, _c in cards)
        return ('<div class="section-header">FITS header</div>'
                f'<div class="meta-table-wrap"><table class="meta-table">{body}</table></div>')
    return ""


def inspector_html(canvas):
    """Details of the object selected in the left sidebar."""
    obj = canvas.selected()
    if obj is None:
        return '<div class="field-static">No selection</div>'
    lo, hi = canvas._display_range(obj.spectrum)
    name = html.escape(obj.name)
    return (
        # the name is editable in place: a change is sent on the "object-rename" channel
        '<div class="field-group"><label>Name</label>'
        f'<input type="text" class="field-input" id="object-rename" value="{name}" '
        f'data-uid="{html.escape(obj.uid)}" spellcheck="false" title="Rename this object">'
        '</div>'
        + _field("Range", f"{lo:.1f} – {hi:.1f} {canvas.xlabel}")
        + _field("Points", f"{len(obj.spectrum):,}")
        + _field("Native units", f"{obj.raw.xunit} / {obj.raw.yunit}")
        + _meta_html(obj.meta)
    )


def footer_html(canvas):
    objs = canvas.objects
    if not objs:
        return '<span class="pill">0 objects</span>'
    lo, hi = canvas.display_bounds()
    rng = f"{lo:.1f} – {hi:.1f} {canvas.xlabel}"
    noun = "object" if len(objs) == 1 else "objects"
    return f'<span class="pill">{len(objs)} {noun}</span><span class="pill">{rng}</span>'


def build():
    # the canvas starts empty — the user opens or synthesizes the first spectrum
    start_dir = picker_start_dir()

    # every hidden bridge / action channel made below is collected here, so the document
    # root column at the end is just `*B.models` instead of a hand-listed enumeration.
    # The context carries it (plus the ToolManager and, later, canvas/notify/run_busy)
    # into the tools, which create and own their bridges themselves.
    B = BridgeSet()
    tools = ToolManager()
    ctx = AppContext(curdoc(), B, tools)

    objects_bridge = B.html("objects")
    # the job list, bottom left; named so a harness can read what the panel shows
    processes_bridge = B.html("processes", name="processes")
    processes_usage_bridge = B.html("processes-usage", name="processes-usage")
    # the job details and recovery modals, filled before they open
    job_details_bridge = B.html("job-details", name="job-details")
    job_modal = B.dialog("job-modal")
    recovery_bridge = B.html("recovery-list", name="recovery-list")
    recovery_modal = B.dialog("recovery-modal")
    inspector_bridge = B.html("inspector")
    footer_bridge = B.html("footer-info")
    # transient info messages (left of footer); named for the same reason
    status_bridge = B.html("footer-status", name="footer-status")

    # show a short message in the footer status slot, then quietly revert to "Ready".
    # the token guards against an earlier message's revert wiping out a newer one.
    # `level` (info/warn/danger/success) drives the semantic colour + the pop/glow flash;
    # each call replaces the inner node, so the CSS animation restarts and re-catches the eye.
    _status_token = itertools.count()
    _status_last = {"n": -1}
    _status_levels = {"info", "warn", "danger", "success"}

    def notify(message, level="warn", revert_ms=4000):
        lvl = level if level in _status_levels else "info"
        n = next(_status_token)
        _status_last["n"] = n
        # trailing comment keeps the pushed string unique, so an identical repeat
        # message still fires the bridge's change event and replays the flash
        status_bridge.value = (f'<span class="status-msg status-{lvl}">'
                               f'{html.escape(str(message))}</span><!--{n}-->')

        def restore():
            if _status_last["n"] == n:
                status_bridge.value = "Ready"

        try:
            ctx.doc.add_timeout_callback(restore, revert_ms)
        except Exception:
            pass

    def remember_recent(path):
        """Add ``path`` to the recent-files list: the Open spectrum picker's history (spectra
        and run folders), which is also where that picker starts. Model atoms are opened
        from their own rooted picker and are not recorded. The list is a convenience, so a config
        the store rejects (unreadable, unsupported or unwritable) is reported in the footer
        and logged, and the open that asked for it goes ahead."""
        try:
            store.add_recent(path)
        except store.ConfigError as exc:
            logging.exception("Could not record %s as a recent file", path)
            notify(f"Recent files not updated: {exc}", level="warn")

    # (the Filter / Normalise / Mask tools create and own their bridges via the context)
    linelist_bridge = B.html("linelist-table")
    linelist_mode_bridge = B.cls("#main", "linelist-mode")  # the line-list bottom panel
    linelist_btn_bridge = B.cls("#line-list-btn", "active")  # its toolbar toggle
    # the panel's search field: its own channel, since a query holds spaces and may hold
    # the ":" the object-actions bridge splits on; NUL-separated from the nonce
    linesearch_actions = B.action("linelist-search")
    linesearch_bridge = B.value("linelist-search")
    run_name_bridge = B.html("run-name")          # name shown on the "Load run" button
    run_ready_bridge = B.cls("#load-run", "ready")  # reveal it inside a run folder
    loading_bridge = B.cls("#footer", "loading")       # footer progress bar while loading
    empty_bridge = B.cls("#main", "empty")  # swap the empty canvas for the inline session loader
    synth_warn_bridge = B.cls("#settings-warning", "show")       # missing-files sign on the cog
    # (the abundance / correction wizards create and own their bridges via the context)

    # which left-sidebar list owns the single selection ("spectrum" | "result" | "corr");
    # the wizards drive it through ctx.focus_kind + the "focus-changed" event, and the
    # Objects list repaints to follow the highlight.
    ctx.focus_kind = "spectrum"
    ctx.subscribe("focus-changed",
                  lambda kind, **_: setattr(objects_bridge, "value",
                                            objects_panel_html(canvas, kind)))

    def synth_on_warn(has_warning):
        toggle(synth_warn_bridge, has_warning)

    inspector = linenav = abundancetool = abundcorrtool = ionbaltool = None  # assigned after the canvas
    atomtool = None
    warm_strengths = None   # ditto: update_all may run before it is defined
    picker = None           # ditto: the spectrum file picker, whose rows spin while loading
    autosave = None         # ditto: the workspace autosave, defined with the session actions
    jobs = ctx.jobs         # every background operation this document admitted

    def run_busy(fn):
        """Defer a UI-only callback. Computational work must use ctx.submit."""
        toggle(loading_bridge, True)
        def work():
            try:
                fn()
            finally:
                toggle(loading_bridge, False)
        curdoc().add_next_tick_callback(work)

    def _range(r):
        """``[start, end]`` of a Bokeh range when it holds an explicit finite interval."""
        lo, hi = r.start, r.end
        if not all(isinstance(v, (int, float)) and np.isfinite(v) for v in (lo, hi)):
            return None
        return [float(lo), float(hi)] if lo < hi else None

    def snapshot():
        """The live workspace as a :class:`session.Workspace` dict."""
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
        focus_kind = ctx.focus_kind
        if focus_kind == "spectrum":
            focus_uid = spectrum_selected
        elif focus_kind == "linemask":
            focus_uid = mask_store.current_uid
        elif focus_kind == "atom":
            focus_uid = atomtool.selected_uid
        else:
            focus_uid = wizards[FOCUS_TOOL[focus_kind]].selected_uid
        return Workspace(
            spectra=specs, spectrum_selected=spectrum_selected,
            masks=MaskSnap(**mask_store.snapshot()),
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
            tools=ToolsSnap(active=tools.active_name(), atom_tab=atomtool.tab,
                            steps={name: w.step for name, w in wizards.items()},
                            jumped={name: w._jumped for name, w in wizards.items()}),
        ).to_dict()

    # the daemon's running and queued jobs and its thread use, once this process holds a
    # daemon connection (a Bokeh server connects when it loads, a harness when it first
    # submits); each change repaints the Processes panel, coalesced per tick
    live_view = {"view": None, "unwatch": None, "repaint": False}

    def live():
        try:
            view = workflow_daemon.live_jobs()
        except (OSError, workflow_daemon.DaemonError):
            log.warning("The workflow daemon cannot be reached", exc_info=True)
            return None
        if view is not None and view is not live_view["view"]:
            if live_view["unwatch"] is not None:
                live_view["unwatch"]()
            live_view["view"], live_view["unwatch"] = view, view.watch(_on_live_event)
        return view

    def _on_live_event(event):
        if ctx.closed or live_view["repaint"]:
            return
        live_view["repaint"] = True

        def apply():
            live_view["repaint"] = False
            repaint_processes()
        ctx.schedule(apply)

    def repaint_processes():
        view = live()
        processes_bridge.value = processes_html(jobs, view, ctx.recovery)
        processes_usage_bridge.value = processes.usage_html(view)

    def update_all(canvas):
        objects_bridge.value = objects_panel_html(canvas, ctx.focus_kind)
        repaint_processes()
        inspector_bridge.value = inspector_html(canvas)
        footer_bridge.value = footer_html(canvas)
        synth_modal.refresh_targets()  # keep the fit target dropdown in sync
        active = bool(jobs.active)
        toggle(loading_bridge, active)
        # show the inline session loader while the canvas holds nothing (and nothing
        # is on its way in; a grid download or listing delivers nothing to the canvas)
        toggle(empty_bridge, not canvas.objects
               and not any(job.delivers for job in jobs.active))
        if inspector is not None:
            inspector.bind(canvas.selected())
        if picker is not None:
            picker.refresh_busy()  # a load settled or was cancelled: stop its row spinning
        ctx.emit("changed")  # EW / Filter / Normalise (and future tools) resync themselves
        if ctx.annot is not None:
            ctx.annot.resync()   # a different object may be a different star (linestrength)
        if warm_strengths is not None:
            warm_strengths()     # ... whose reference list may still need reading
        if linenav is not None and linenav.active:
            linenav.render()  # keep the line list scoped to the active object's range
            # (only while the panel is open; it is rebuilt on the next open otherwise)
        if abundancetool is not None:
            abundancetool.refresh()  # keep the abundance modal's reference selector in sync
        if ionbaltool is not None:
            ionbaltool.refresh()     # ditto for the balance modal's target selector
        if autosave is not None:
            autosave()               # the workspace changed: keep the autosave slot current

    canvas = SpectrumCanvas(on_change=update_all)  # empty until a spectrum is opened
    ctx.canvas = canvas
    ctx.notify = notify
    ctx.run_busy = run_busy
    ctx.busy = lambda on: toggle(loading_bridge, bool(on))

    def edit_fit():
        """Edit-tab "Fit to target": broaden the selected (model) spectrum to best match
        another object over the current view — pure post-processing, no re-synthesis."""
        obj = canvas.selected()
        if obj is None:
            return
        target = canvas._find(inspector.target)
        if target is None or target is obj:
            notify("Pick another spectrum as the fit target.")
            return
        params = inspector.armed()
        if not params:
            notify("Switch on at least one parameter to fit.")
            return
        r = canvas.main.x_range
        if r.start is None or r.end is None:
            notify("Set a view range to fit over.")
            return
        lo, hi = to_range(r.start, r.end, canvas.xunit, target.spectrum.xunit)

        source, observed, initial = obj.spectrum, target.spectrum, dict(obj.conv)

        def run(token):
            token.check()
            tgt = observed.crop(xmin=lo, xmax=hi)
            best, _ = source.fit(tgt, fit=params, p0=initial)
            token.check()
            return best

        def done(best):
            if canvas._find(obj.uid) is not obj or obj.conv != initial:
                return  # removed or edited while the fit was running
            obj.conv = {k: v for k, v in {**initial, **best}.items() if v}
            canvas._refresh_object(obj)
            if canvas.selected() is obj:
                inspector.bind(obj)
            update_all(canvas)
            notify("Broadening fit · "
                   + ", ".join(f"{k}={v:g}" for k, v in best.items()), level="success")

        load_async("broadening-fit", run, done, executor=_SYNTH_EXEC)

    def edit_mask():
        """Edit-tab "Select full range": the Select range tool on the selected spectrum's extent."""
        obj = canvas.selected()
        if obj is None:
            return
        x, _ = obj.spectrum.values
        lo, hi = to_range(float(x.min()), float(x.max()), obj.spectrum.xunit, u.AA)
        range_tool.select_range(lo, hi)

    inspector = Inspector(canvas, on_fit=edit_fit, on_mask=edit_mask, notify=notify)
    # the Select range tool: the one range gesture, whose actions reach the line-mask view and the
    # synthesis modal built further down (the lambdas resolve them at click time)
    range_tool = RangeTool(
        ctx,
        on_new_mask=lambda lo, hi: masktool.create_from_window(lo, hi),
        on_synth=lambda lo, hi: synth_modal.open_range(lo, hi, "synth"),
        on_fit=lambda lo, hi: synth_modal.open_range(lo, hi, "fit"))


    prefs = load_prefs()
    master = str(_resolve_path(prefs, prefs["linelist"]))
    # the editable reference-star catalogue (Preferences → Reference stars) is installed
    # into tumult.refstars here, before anything lists it; its edits repaint every dropdown
    # that does (on_refstars_change, defined once those exist)
    refstars_page = RefStarsPage(on_change=lambda: on_refstars_change())
    strengths = LineStrengths(master_path=master, linelist=master,
                              refstar=store.get_pref("strength_refstar", refstars.AUTO_KEY))
    ctx.strengths = strengths

    linenav = LineNavigator(
        canvas, strengths, on_table=lambda h: setattr(linelist_bridge, "value", h),
        on_open=lambda on: (toggle(linelist_mode_bridge, on), toggle(linelist_btn_bridge, on)),
        on_query_cleared=lambda: setattr(linesearch_bridge, "value", f"\0{time.time()}"))
    linesearch_actions.on_change(
        "value", lambda attr, old, new: linenav.set_query(new.split("\0")[0]))

    # renaming from the Info tab: "<uid>\0<name>\0<nonce>" (names may hold spaces/colons)
    rename_actions = B.action("object-rename")

    def on_rename(attr, old, new):
        if not new:
            return
        uid, name = new.split("\0")[:2]
        canvas.rename_object(uid, name)

    rename_actions.on_change("value", on_rename)

    def on_line_setting(key, value):
        if key == "refstar":
            # a different star is a different ranking, so the overlay and the table both have
            # to be redrawn from it. Repainting *now* would only show the proxy — the new
            # band is not computed yet — so when there is something to compute, leave the
            # previous star's names up and let warm_strengths repaint once the real ranking
            # is there. Flashing the proxy's TiO and C2 in between is worse than a short wait.
            store.set_pref("strength_refstar", strengths.set_refstar(value))
            if not warm_strengths():
                ctx.annot.resync()
                linenav.render()
        else:
            linenav.set_jump_distance(value)

    # the line-list panel's own controls, rendered inline in its header
    linejump_form = Form("linesettings", [
        # which atmosphere the strengths — and so the ranking, and the order of this table —
        # are computed on. Named classes rather than a Teff: it is what a user of an
        # unlabelled observation actually knows about their star (see refstars)
        Field("refstar", "Reference star", kind="select", value=strengths.refstar,
              options=refstars.options(), inline=True),
        Field("dist", "Distance from center [Å]", kind="number",
              value=linenav.jump_distance, min=0.01, inline=True),
    ], on_change=on_line_setting)
    # line-position indicators over the spectrum, shared by whichever view wants them (the
    # LineMasks view shows its set's element; other tools attach the same way). Unowned, it
    # draws the strongest lines of the window — as many as the frame is wide enough to name.
    ctx.annot = LineAnnotations(ctx, strengths)
    # A pan or zoom can leave the computed band; ask for the new one (a no-op when it is
    # already cached, so this is free on the overwhelming majority of frames). Continuous
    # scrolling crosses band after band, and each crossing would otherwise start a run that
    # the next crossing makes irrelevant, so the request waits out the scroll: every frame
    # replaces the pending one, and only the band the view comes to rest in is computed.
    view_delay = Trailing(STRENGTH_DELAY_MS)
    canvas.add_view_listener(lambda: view_delay("strengths", warm_strengths))

    strength_pending = set()

    def warm_strengths():
        """Compute the strengths for the visible band, off the document thread.

        A band is an EOS solve plus an opacity pass — a fraction of a second for the usual
        master list, seconds for a wide view or a folder of per-species lists — so it must
        not run on a repaint. Until it lands the overlay keeps whatever it last drew, and
        this repaints once the real thing is there.

        Returns whether it scheduled anything, so a caller that changed the ranking can
        repaint immediately when the band was already cached and there is nothing to wait for.
        """
        if not canvas.objects:
            return False
        # While a tool owns the overlay it draws its own line frame, not the global
        # ranking — nothing to warm unless the line-list panel is also open and reading
        # that ranking directly. Otherwise a tool with no dump of its own (an NLTE
        # correction run without a linelist, say) attaches windows that rarely match an
        # already-cached band, and every click would queue a real EOS-and-opacity job.
        if ctx.annot.owner is not None and not linenav.active:
            return False
        obj = canvas.selected()
        lo, hi = ctx.annot._window_aa()
        if not strengths.pending(obj, lo, hi):
            return False

        # Snapshot the reference and band on the document thread. Do not read
        # canvas objects or form values from the HOP worker.
        star, linelist, band = strengths.request_for(obj, lo, hi)
        request = (star.key, linelist, band)
        if request in strength_pending:
            return False
        strength_pending.add(request)

        def run(token):
            token.check()
            return linestrength.compute_band(star, linelist, *band)

        def done(frame):
            ctx.annot.resync()
            if linenav.active:
                linenav.render()
            if frame is None:
                notify("Line strengths unavailable — ranking is approximate", level="warn")

        # Schedule after this callback, avoiding recursive update_all during submit.
        def submit_warm():
            future = load_async("line-strengths", run, done, executor=_HOP_EXEC)
            if future is None:          # refused by the pool: the row says why
                strength_pending.discard(request)
                return
            future.add_done_callback(
                lambda _: ctx.schedule(lambda: strength_pending.discard(request)))
        ctx.schedule(submit_warm)
        return True

    filtertool = FilterTool(ctx)
    normtool = NormaliseTool(ctx)
    # the shared LineMask set (imported windows for fitting / local renorm), owned centrally so
    # the canvas and — later — the abundance / synth-fit workflows all read one source of truth.
    _linemask_dir = (Path(str(store.get_pref("dispatch_cwd") or HERE.parent)).expanduser()
                     / "input" / "linemasks").resolve()
    mask_store = MaskStore(ctx)
    ctx.masks = mask_store
    masktool = MaskTool(
        ctx, mask_store,
        input_dir=_linemask_dir if _linemask_dir.is_dir() else Path.home(),
        # the Auto-generate source runs the wizards' diagnostic line selection over a loaded
        # spectrum (no synthesis), so it needs the canvas objects and a line list
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master, profiles=PROFILES,
        # "+ New mask" draws its window with the Select range tool's gesture
        place_window=range_tool.begin_add,
        # the mask view's analysis actions: hand the open set to the abundance / balance
        # tools and the synthesis modal, which are built further down (the lambdas resolve
        # them at click time)
        on_fit_balance=lambda mset: ionbaltool.open_for(mset),
        on_measure_abund=lambda mset: abundancetool.open_for(mset),
        on_synth=lambda mset: synth_modal.open_masks(mset, "synth"),
        on_fit=lambda mset: synth_modal.open_masks(mset, "fit"))
    settings = build_settings(canvas)

    def add_synth(spec, name, meta=None, conv=None):
        color, dash = canvas.next_run_style()  # distinct colour, picked at add time
        # data embedded on save (source path None); a fit also passes its fitted broadening
        # as `conv`, applied live so the Edit tab shows it and the stored spectrum stays raw.
        o = canvas.add_object(spec, name=name, color=color, dash=dash,
                              source=("synth", None), meta=meta)
        if conv:
            o.conv = dict(conv)
            canvas._refresh_object(o)  # update_all (in _settle) then rebinds the Edit tab
        synth_modal.refresh_cache()  # the new run is now in the cache listing

    # every synthesis goes through ctx.submit_workflow and every fit through
    # ctx.submit_task (durable, detached); both are installed below, before any
    # submission can run
    synth_modal = SynthModal(ctx, add=add_synth, on_warn=synth_on_warn)
    ctx.synth_modal = synth_modal      # local harnesses drive the modal through it

    def marcs_for_strengths():
        """Name the registered MARCS grid the line-strength reference atmospheres are read
        from (the first one in registry order), or none."""
        grid = next((g for g in synth_modal.grids.grids
                     if str(g.get("format") or "").lower() == "marcs"
                     and (g.get("folder") or "").strip()), None)
        linestrength.set_marcs_grid(grid["folder"] if grid else None,
                                    grid.get("pattern") if grid else None)

    marcs_for_strengths()

    # the Abundance tool owns its synthesis setup (its own Setup sub-tab), sharing only the
    # app-global grid registry and the synthesis modal's NLTE model-atom selection; it inherits
    # the modal's configuration on first open and is independent thereafter.
    # its fields live in the abundance modal, so every edit restates the setup in the object
    # panel, where the Run synthesis button stands
    _abund_note = lambda *_: (abundancetool is not None and abundancetool.push_setup_note())
    from atompicker import AtomPicker
    abund_atoms = AtomPicker(Path(load_prefs()["dispatch_cwd"]) / "input" / "atoms", ns="abund-atoms", section="abund_run",
        picker_sel="#abund-atom-picker", section_sel="#abund-atoms",
        actions_name="abund-atom-actions", results_id="abund-atom-results",
        layout_name="abund_atom_picker", del_class="abund-atom-del")
    abund_setup = SynthSetup(
        "abund", "abund_run", synth_modal.grids,
        atoms_provider=lambda: list(abund_atoms.selected), notify=notify,
        on_run_change=_abund_note, on_base_change=_abund_note,
        on_mode_change=lambda on: (abund_atoms.show(on), _abund_note()), on_atmos_change=_abund_note)
    abundancetool = AbundanceTool(
        ctx,
        setup=abund_setup,
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        # follow the canvas's Filtered switch: when filtering is on, measure/fit the filtered
        # spectrum (so EW hatches + line fits match what is shown), else the raw observation
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master)

    # the abundance-correction tool: its own synthesis setup (target/NLTE atmosphere) in its
    # modal, sharing only the app-global grid registry (one grid flagged as the 1D LTE
    # reference). The model atom is the atom object the modal was opened from.
    corr_setup = SynthSetup(
        "corr", "abundcorr_run", synth_modal.grids,
        atoms_provider=lambda: ([abundcorrtool._pending_atom.path]
                                if abundcorrtool._pending_atom is not None else []),
        notify=notify)

    def add_corr_spectrum(spec, name):
        color, dash = canvas.next_run_style()
        canvas.add_object(spec, name=name, color=color, dash=dash, source=("synth", None))

    abundcorrtool = AbundCorrTool(
        ctx,
        setup=corr_setup,
        grids=synth_modal.grids,
        add_spectrum=add_corr_spectrum)

    # the model-atom objects: Open → Model atom loads one; its view is the two
    # tables, its panel the entry point of the correction modal, and its correction results
    # are listed under it
    atomtool = AtomTool(ctx, corr=abundcorrtool)
    abundcorrtool.on_results = atomtool.render_objects
    abundcorrtool._known_atom_paths = lambda: [a.path for a in atomtool.atoms]
    _atoms_dir = Path(str(store.get_pref("dispatch_cwd") or ".")).expanduser() / "input" / "atoms"

    def on_atom_select(filepath):
        from tumult.atom import AtomClass
        name = Path(filepath).name
        load_async(name, lambda token=None: AtomClass(str(filepath)),
                   lambda atom: atomtool.add(filepath, atom), path=filepath)

    # The ionisation-balance modal owns its grid selection and always interpolates it.  Only
    # the app-global registered-grid list and synthesis preferences are shared.
    ionbaltool = IonBalanceTool(
        ctx,
        grids=synth_modal.grids,
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master)

    def on_grids_change():
        """The grid registry changed: every setup re-lists its grid, and the strength
        ranking follows the MARCS grid."""
        abund_setup.refresh_grids()
        corr_setup.refresh_grids()
        marcs_for_strengths()
        warm_strengths()

    synth_modal.on_grids_change = on_grids_change

    def on_refstars_change():
        """The reference-star catalogue was edited: repaint every dropdown listing it, and
        rerank by the edited star if it is the one in use."""
        for setup in (synth_modal.setup, abund_setup, corr_setup):
            setup.refresh_refstars()
        ionbaltool.refresh_refstars()
        field = linejump_form.fields["refstar"]
        field.options = refstars.options()
        field.value = strengths.set_refstar(strengths.refstar)
        store.set_pref("strength_refstar", field.value)
        linejump_form.render()
        if not warm_strengths():
            ctx.annot.resync()
            linenav.render()

    # single-selection coordination between the Objects and the wizards' result lists
    def focus_spectrum(uid):
        if tools.active_name() in ("abundance", "abundcorr", "ionbal", "linemask", "atom"):
            tools.deactivate_active()  # selecting a spectrum leaves the wizard / mask view
        abundancetool.blur()          # leave the centre abund view, drop the result highlight
        abundcorrtool.blur()          # ditto for the correction view
        ionbaltool.blur()             # ditto for the measurement view
        canvas.select_object(uid)     # -> on_change -> update_all repaints with spectrum focus

    # a measurement's panel names the set it came from (and flags it once that set is edited)
    ctx.subscribe("masks-changed",
                  lambda **_: ionbaltool.active and ionbaltool.push_panel())

    # the abundance / parameter wizards' step-1 EWs track the reference's filter state
    ctx.subscribe("filter-state", lambda **_: abundancetool.on_filter_changed())
    ctx.subscribe("filter-state", lambda **_: ionbaltool.on_filter_changed())

    # the mutually-exclusive sidebar tools: at most one owns the right panel at a time.
    # filter's enter is wrapped in run_busy so the progress bar paints during the FFT/fit.
    tools.register("range", range_tool)
    tools.register("normalise", normtool)
    tools.register("linemask", masktool)
    tools.register("filter", filtertool, runner=run_busy)
    tools.register("abundance", abundancetool)
    tools.register("abundcorr", abundcorrtool)
    tools.register("ionbal", ionbaltool)
    tools.register("atom", atomtool)

    # --- session save / restore -------------------------------------------
    # the three result lists, by the name of their tool (the session's result sections)
    wizards = {"abundance": abundancetool, "abundcorr": abundcorrtool, "ionbal": ionbaltool}
    assert tuple(wizards) == RESULT_TOOLS

    def _restore_range(rng, vals):
        if not vals:
            return
        lo, hi = vals
        if rng.bounds:
            lo, hi = max(lo, min(rng.bounds)), min(hi, max(rng.bounds))
        if lo < hi:
            rng.start, rng.end = lo, hi

    def prepare_workspace(snap):
        """Validate a slot and build its objects off the live workspace
        (:func:`session.prepare`); :class:`SessionError` when the slot is unusable."""
        from tumult.atom import AtomClass
        ws = Workspace.from_dict(snap)
        return session.prepare(ws, load_spectrum=_load_obj,
                               load_atom=lambda path: AtomClass(str(path)), tools=wizards)

    def _capture_live():
        """The live workspace as a prepared one (no file reads): what a failed commit
        re-applies."""
        ws = Workspace.from_dict(snapshot())
        by_uid = {o.uid: o for o in canvas.objects}
        return PreparedWorkspace(
            spectra=[PreparedSpectrum(sp, by_uid[sp.uid].raw, by_uid[sp.uid].meta)
                     for sp in ws.spectra],
            spectrum_selected=ws.spectrum_selected,
            mask_sets=list(mask_store.sets), mask_current=mask_store.current_uid,
            mask_active=mask_store.active,
            atoms=[PreparedAtom(AtomSnap(uid=a.uid, path=a.path), a.atom)
                   for a in atomtool.atoms],
            atom_selected=atomtool.selected_uid,
            results={name: (list(w.results), w.selected_uid) for name, w in wizards.items()},
            view=ws.view, filter=ws.filter, focus=ws.focus, tools=ws.tools)

    def _place_spectrum(ps):
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

    def _apply_workspace(p):
        """Replace every collection with the prepared workspace ``p`` and restore its
        view, selections and focus. Every collection is replaced, empty ones included;
        spectra, masks and atoms land before the results that refer to them."""
        tools.deactivate_active()           # whatever view was open, the canvas returns
        for w in wizards.values():
            w.blur()
        ctx.focus_kind = "spectrum"
        filtertool.clear()
        for o in list(canvas.objects):
            canvas.remove_object(o.uid)
        for ps in p.spectra:
            _place_spectrum(ps)
        canvas.selected_uid = p.spectrum_selected     # None when none was, or it is gone
        mask_store.restore(p.mask_sets, p.mask_current, p.mask_active)
        atomtool.atoms = [AtomObject(path=a.snap.path, atom=a.atom, uid=a.snap.uid)
                          for a in p.atoms]
        atomtool.selected_uid = p.atom_selected
        atomtool.set_tab(p.tools.atom_tab)
        for name, w in wizards.items():
            results, selected = p.results[name]
            w.restore(results, selected)
            if name in p.tools.steps:
                w.set_step(p.tools.steps[name])
        atomtool.render_objects()
        v = p.view
        if v.xunit:
            canvas.set_xunit(v.xunit)
        if v.yunit:
            canvas.set_yunit(v.yunit)
        settings.set_values(xunit=canvas.xlabel, yunit=canvas.ylabel)  # sync the prefs form
        f = p.filter
        if canvas.objects and f.method == "lowpass" and f.res is not None:
            filtertool.apply_saved_lowpass(f.res)
            if not f.active:
                filtertool.set_visible(False)
        elif canvas.objects and f.fmax is not None and f.floor is not None:
            filtertool.apply_saved(f.fmax, f.floor)
            if not f.active:
                filtertool.set_visible(False)
        # the focus, through the same routing a click takes, so the panels it opens are
        # the ones the user had; a cleared focus leaves the plain canvas
        focus = p.focus
        if focus.kind == "spectrum":
            if focus.uid is not None:
                focus_spectrum(focus.uid)
        elif focus.kind == "linemask":
            masktool.select_set(focus.uid)
        elif focus.kind == "atom":
            atomtool.select(focus.uid)
        elif focus.uid is not None:
            wizards[FOCUS_TOOL[focus.kind]].select_result(focus.uid)
        for name, w in wizards.items():
            jumped = p.tools.jumped.get(name)
            res = w.selected()
            if jumped is None or res is None or jumped >= len(res.lines):
                continue
            if w.active and w.is_focused():
                w._jump_line(str(jumped))      # the accent, hatch and CoG of that line
            else:
                w._jumped = jumped             # painted when its view is next opened
        # a panel open beside the canvas (a wizard's steps, the mask or atom view) is
        # reopened; the range, filter and normalise gestures are transient and stay closed
        active = p.tools.active
        if active in ("abundance", "abundcorr", "ionbal", "linemask", "atom") \
                and tools.active_name() != active:
            tools.activate(active)
        # the explicit ranges last, after everything above that moves the view
        _restore_range(canvas.nav.x_range, v.nav_x)
        _restore_range(canvas.main.x_range, v.main_x)
        _restore_range(canvas.main.y_range, v.main_y)
        _restore_range(canvas.nav.y_range, v.nav_y)
        update_all(canvas)

    def restore_workspace(snap, name="save", workflows=()):
        """Load a slot as one document operation. The slot is validated and its objects
        built first (a broken slot changes nothing); the commit then marks a workspace
        boundary (new generation, old result subscriptions forgotten, this document's
        non-detached jobs cancelled) and replaces every collection. If applying fails,
        the previous workspace is put back. The slot's ``workflows`` references then
        replace this document's (see ``_reattach``). Returns the
        :class:`session.LoadReport`, or ``None`` when nothing was loaded."""
        try:
            prepared = prepare_workspace(snap)
            live = _capture_live()
        except Exception as exc:      # SessionError names the field; anything else is logged
            logging.exception("The saved workspace %r cannot be loaded", name)
            notify(f"Could not load “{name}”: {exc}", level="danger", revert_ms=10000)
            return None
        ctx.advance_generation()
        cancelled = [job.name for job in jobs.active if not job.detached]
        for job in list(jobs.active):
            if not job.detached:
                jobs.request_cancel(job.pid)
        try:
            _apply_workspace(prepared)
        except Exception as exc:
            logging.exception("Applying the saved workspace %r failed; restoring the "
                              "previous workspace", name)
            try:
                _apply_workspace(live)
            except Exception:
                logging.exception("The previous workspace could not be restored either")
                notify(f"Loading “{name}” failed: {exc}; the previous workspace could not "
                       "be restored", level="danger", revert_ms=10000)
                return None
            notify(f"Loading “{name}” failed: {exc}; the previous workspace was restored",
                   level="danger", revert_ms=10000)
            return None
        outcomes, gone = _reattach(workflows)
        report = session.LoadReport(name, prepared.counts(), prepared.skipped + gone,
                                    prepared.unresolved, cancelled, outcomes)
        update_all(canvas)
        for line in report.details():
            log.warning("%s", line)
        notify(report.message(), level="warn" if report.partial else "success",
               revert_ms=12000 if report.partial else 4000)
        return report

    ctx.restore_workspace = restore_workspace   # local harnesses load slots through it

    # --- background work: a Processes row now, the outcome delivered on the document ----
    def _admit(name, *, heavy=False, detached=False, path=None, delivers=True):
        job = Job(_nonce(), name, synth.CancelToken(), ctx.generation,
                  heavy=heavy, detached=detached, path=str(path) if path else None,
                  delivers=delivers)
        jobs.add(job)
        # the row appears once the job has lasted SHOW_AFTER_S (see processes_html)
        try:
            ctx.doc.add_timeout_callback(lambda: None if job.terminal or ctx.closed
                                         else repaint_processes(),
                                         int(processes.SHOW_AFTER_S * 1000) + 20)
        except Exception:
            pass                     # the document went away
        return job

    def _admission_failed(job, exc):
        """The pool refused the job: the same terminal path as a failure in flight, so
        the spinner clears and the row keeps the reason."""
        log.exception("%s (job %d) could not be admitted", job.name, job.pid)
        jobs.settle(job.pid, FAILED, error=f"{type(exc).__name__}: {exc}")
        notify(f"{job.name} could not start: {exc} — see Processes", level="danger")
        update_all(canvas)

    def load_async(name, loader, on_done, executor=_EXEC, path=None, detached=False,
                   delivers=True):
        """Run ``loader(token)`` on ``executor`` and hand its value to ``on_done`` on the
        document thread. Returns the pool future, or ``None`` when the pool refused it.
        ``delivers`` False marks work whose outcome stays off the canvas (a grid
        download), so the empty state does not wait for it."""
        job = _admit(name, heavy=executor is _SYNTH_EXEC, detached=detached, path=path,
                     delivers=delivers)
        token = job.token                # the loader threads it into synth.* for hard cancel
        update_all(canvas)               # show the queued row immediately

        def work():                      # runs in a worker thread
            if token.cancelled:          # cancelled before the worker even started
                _settle(job, None, None, None)
                return
            if job.state == QUEUED:
                job.state, job.started = RUNNING, time.time()
                ctx.schedule(lambda: update_all(canvas))
            try:
                result = profiled(name, lambda: loader(token))
            except synth.CancelledError:
                _settle(job, None, None, None)  # token.cancelled -> dropped
            except Exception as exc:
                # the full traceback belongs to the server log whether or not the
                # document that asked still exists; the row gets the concise message
                log.exception("%s (job %d) failed", name, job.pid)
                _settle(job, None, None, exc)
            else:
                _settle(job, result, on_done, None)

        try:
            job.future = executor.submit(work)
        except Exception as exc:
            _admission_failed(job, exc)
            return None
        return job.future

    # Document-owned work on the CPU-budgeted pool (a line measurement, a broadening fit,
    # a normalisation): cancelled when the session closes. Dispatch calculations do not
    # come through here; they are durable workflows (submit_workflow, submit_task).
    ctx.submit = lambda name, loader, on_done, detached=False: load_async(
        name, loader, on_done, executor=_SYNTH_EXEC, detached=detached)
    # network transfers (a grid download) wait on I/O, not the CPU budget: the IO pool
    ctx.submit_io = lambda name, loader, on_done: load_async(name, loader, on_done,
                                                            delivers=False)
    ctx.execution = _get_runtime()

    def _import_result(handle):
        """The delivery of a workflow this document did not submit: a synthesis adds
        every spectrum of its result under the name and annotations recorded at
        submission; a tool calculation goes to the importer its tool registered."""
        kind = handle.kind.name
        if kind != "synthesis":
            importer = ctx.task_importers.get(kind)
            if importer is None:
                raise ValueError(f"no tool imports a {kind} result")
            return lambda result: importer(result, handle)

        def deliver(result):
            for index, spec in enumerate(result.spectra, start=1):
                suffix = f"·{index}" if len(result.spectra) > 1 else ""
                add_synth(spec, handle.label + suffix, meta=handle.meta)
        return deliver

    def _import_progress(handle):
        """The progress delivery of a tool record this document did not submit, or
        ``None`` for a kind without progress or without a registered tool."""
        kind = handle.kind
        follow = ctx.task_progress.get(kind.name) if kind.progress else None
        return None if follow is None else (lambda snapshot: follow(snapshot, handle))

    def _subscribe(job, handle, on_done, on_error=None):
        """Deliver ``handle``'s outcome to ``job`` once it settles; the delivery marks
        the workflow imported so a later load does not offer its result again."""
        def deliver(result):
            on_done(result)
            ctx.workflow_refs.mark_imported(handle.id)

        def settled(handle):             # worker thread (or the canceller's), never waits
            try:
                result = handle.result()
            except synth.CancelledError:
                _settle(job, None, None, None)
            except Exception as exc:
                _settle(job, None, None, exc, on_error)
            else:
                _settle(job, result, deliver, None)

        handle.add_done_callback(settled)

    PROGRESS_POLL_MS = 2000        # Slurm records only; the daemon reports its snapshots

    def _follow_progress(job, handle, on_progress):
        """While ``job`` runs, hand every new progress snapshot of its record to
        ``on_progress`` on the document thread. The daemon reports each snapshot's stamp
        (a Slurm record is polled instead); the snapshot is decoded off the document
        thread, since it can carry a spectrum, one read at a time, and a stamp that
        arrives during a read is read next. Ends with the job, the document or the
        workspace generation."""
        if on_progress is None or not handle.kind.progress:
            return
        state = {"stamp": None, "busy": False, "stop": None}

        def stop():
            stopper, state["stop"] = state["stop"], None
            if stopper is not None:
                try:
                    stopper()
                except Exception:
                    pass                 # already removed, or the document is gone

        def tick():
            if job.terminal or ctx.closed or job.generation is not ctx.generation:
                stop()
                return
            if state["busy"]:
                return
            try:
                stamp = handle.progress_stamp()
            except OSError:
                return
            if stamp is None or stamp == state["stamp"]:
                return
            state["busy"] = True

            def read():                  # IO pool
                try:
                    snapshot = handle.progress_result()
                except Exception:
                    log.exception("%s (job %d): progress snapshot unreadable", job.name, job.pid)
                    snapshot = None

                def deliver():
                    state["busy"], state["stamp"] = False, stamp
                    if snapshot is None or job.terminal or job.generation is not ctx.generation:
                        return
                    try:
                        on_progress(snapshot)
                    except Exception:
                        log.exception("%s (job %d): progress not applied", job.name, job.pid)
                    tick()               # a newer snapshot written during the read
                ctx.schedule(deliver)
            _EXEC.submit(read)

        if handle.detached == "slurm":
            callback = ctx.doc.add_periodic_callback(tick, PROGRESS_POLL_MS)
            state["stop"] = lambda: ctx.doc.remove_periodic_callback(callback)
        else:
            state["stop"] = handle.on_progress(lambda stamp: ctx.schedule(tick))
            ctx.schedule(tick)           # a snapshot written before this document followed

    def submit_workflow(name, request, on_done, meta=None):
        """Queue one durable synthesis workflow and deliver its result to this document
        once its handle settles. Submission happens here, once; nothing waits on the
        workflow from another pool thread, and the document thread never blocks on the
        result. The workflow is detached: closing the tab leaves it running, and only its
        Processes row cancels it, through the handle. ``name`` and ``meta`` go into the
        record, so a later document lists and imports it as it was submitted, and the
        workflow is referenced beside the workspace from now on. Returns the handle, or
        ``None`` when submission failed (the row keeps the reason)."""
        job = _admit(name, heavy=True, detached=True)
        # the row's cancel token forwards to the handle, so a settle can tell a cancelled
        # workflow's late outcome from a wanted one the same way it does for loads
        try:
            handle = submit_synthesis(request, label=name, meta=meta)
        except Exception as exc:
            _admission_failed(job, exc)
            return None
        job.handle = handle
        job.token.add_callback(handle.cancel)
        ctx.workflow_refs.attach(handle, name)
        update_all(canvas)
        _subscribe(job, handle, on_done)
        return handle

    ctx.submit_workflow = submit_workflow

    def submit_task(name, kind, args, on_done, *, prefs, meta=None, on_progress=None,
                    on_error=None):
        """Queue one tool calculation (:mod:`tumult.synthesis.tasks`) as a durable
        workflow, as :func:`submit_workflow` does a synthesis: detached, referenced
        beside the workspace, its result delivered to ``on_done`` once, a failure to
        ``on_error`` after the row has recorded it, and every progress snapshot of the
        record to ``on_progress`` while it runs. ``prefs`` are the synthesis preferences
        the calculation was configured with (the driver and the record root follow
        them). Returns the job, or ``None`` when submission failed (the row keeps the
        reason)."""
        job = _admit(name, heavy=True, detached=True)
        try:
            handle = submit_task_record(kind, args, prefs=prefs, label=name, meta=meta)
        except Exception as exc:
            _admission_failed(job, exc)
            return None
        job.handle = handle
        job.token.add_callback(handle.cancel)
        ctx.workflow_refs.attach(handle, name)
        update_all(canvas)
        _subscribe(job, handle, on_done, on_error)
        _follow_progress(job, handle, on_progress)
        return job

    ctx.submit_task = submit_task

    # --- recovery of durable workflows after a reconnect ------------------------------
    def _offer(handle, referenced=False):
        """List ``handle`` for recovery unless this document already follows it or has
        imported its result. Returns the entry, or ``None``."""
        if jobs.find_workflow(handle.id) is not None or ctx.workflow_refs.imported(handle.id):
            return None
        try:
            entry = Recoverable.from_handle(handle, referenced)
        except (OSError, ValueError, KeyError) as exc:
            log.warning("Workflow record %s cannot be read: %s", handle.folder, exc)
            return None
        return ctx.recovery.offer(entry)

    def discover_workflows(limit=20):
        """Ask the daemon for the newest settled records under the configured roots and
        offer the ones this document has not imported in the recovery modal, which opens
        on them. Running and queued workflows are rows of the panel already. Returns the
        ids offered."""
        try:
            workflow_daemon.client()
            view = live()
            reply = view.refresh(workflow_roots(), limit) if view is not None else {"jobs": []}
        except (OSError, workflow_daemon.DaemonError) as exc:
            notify(f"The workflow daemon cannot be reached: {exc}", level="danger")
            return []
        offered = []
        for summary in reply["jobs"]:
            if summary["type"] != "workflow" or summary["state"] in (QUEUED, RUNNING):
                continue
            handle = submission.open_workflow(summary["folder"])
            entry = _offer(handle, referenced=handle.id in ctx.workflow_refs)
            if entry is not None:
                offered.append(entry.id)
        show_recovery()
        return offered

    def show_recovery():
        recovery_bridge.value = processes.recovery_html(ctx.recovery)
        repaint_processes()
        pulse(recovery_modal, "open")

    def _live_summary(job_id):
        view = live()
        return view.get(job_id) if view is not None else None

    def _external_run(run_id):
        summary = _live_summary(run_id)
        if summary is None or summary["type"] != "external":
            return None
        return ExternalRun.from_summary(summary)

    def cancel_run(run_id):
        """End an external run through the daemon, which signals it only when the process
        table confirms the recorded pid still runs that input (or scancels it)."""
        entry = _external_run(run_id)
        if entry is None:
            return
        try:
            state = workflow_daemon.client().request("cancel", job=run_id)["state"]
        except (OSError, workflow_daemon.DaemonError) as exc:
            notify(f"{entry.name} could not be ended: {exc}", level="danger")
            return
        if state == CANCELLED:
            notify(f"Ended {entry.name}", level="info")
            pulse(job_modal, "close")
        else:
            notify(f"{entry.name} could not be signalled from here; it is {state}", level="warn")
        repaint_processes()

    def cancel_workflow(workflow_id):
        """Cancel a daemon workflow this document does not follow (another tab's, a
        script's) from its row."""
        summary = _live_summary(workflow_id)
        folder = summary["folder"] if summary else None
        if folder is None:
            entry = ctx.recovery.get(workflow_id)
            folder = entry.handle.folder if entry is not None else None
        if folder is None:
            return
        handle = submission.open_workflow(folder)
        try:
            told = handle.cancel()
        except (OSError, workflow_daemon.DaemonError) as exc:
            notify(f"{handle.label} could not be cancelled: {exc}", level="danger")
            return
        notify(f"Cancelling {handle.label}…" if told else f"{handle.label} is {handle.status}",
               level="info" if told else "warn")
        pulse(job_modal, "close")
        repaint_processes()

    def use_run_departures(run_id):
        """Hand a departure solve to the synthesis modal's Saved departures field."""
        entry = _external_run(run_id)
        if entry is None:
            return
        synth_modal.setup.mode_form.set_values(mode="NLTE", departure_field=entry.folder)
        synth_modal.setup._save_mode("mode", "NLTE")
        synth_modal.setup._save_mode("departure_field", entry.folder)
        notify(f"Synthesis modal set to reuse departures from {entry.name}", level="info")

    def _details_target(uid):
        """What a details button names: ``("job", Job)``, ``("workflow", handle,
        summary, recoverable)`` or ``("run", ExternalRun)``; ``None`` if it is gone."""
        if uid.startswith("job-"):
            try:
                job = jobs.get(int(uid[4:]))
            except ValueError:
                return None
            return None if job is None else ("job", job)
        summary = _live_summary(uid)
        if summary is not None and summary["type"] == "external":
            return ("run", ExternalRun.from_summary(summary))
        recoverable = ctx.recovery.get(uid)
        if summary is not None:
            return ("workflow", submission.open_workflow(summary["folder"]), summary, None)
        if recoverable is not None:
            return ("workflow", recoverable.handle, None, recoverable)
        ref = ctx.workflow_refs.get(uid)
        if ref is not None:
            return ("workflow", submission.open_workflow(ref.folder), None, None)
        return None

    def show_details(uid):
        """Fill and open the job details modal for a row's details button."""
        target = _details_target(uid)
        if target is None:
            notify("That job is no longer listed", level="warn")
            repaint_processes()
            return
        try:
            if target[0] == "run":
                body = processes.run_details_html(target[1])
            elif target[0] == "job" and target[1].handle is None:
                body = processes.job_details_html(target[1])
            elif target[0] == "job":
                job = target[1]
                body = processes.workflow_details_html(
                    job.handle, _live_summary(job.workflow_id), job=job)
            else:
                _, handle, summary, recoverable = target
                follower = jobs.find_workflow(handle.id)
                body = processes.workflow_details_html(handle, summary, job=follower,
                                                       recoverable=recoverable)
        except (OSError, ValueError, KeyError) as exc:
            notify(f"The details cannot be read: {exc}", level="danger")
            return
        job_details_bridge.value = body
        pulse(job_modal, "open")

    def open_folder(uid):
        """Show a job's folder in the file manager of the machine the server runs on."""
        import subprocess
        target = _details_target(uid)
        if target is None:
            return
        if target[0] == "run":
            folder = target[1].folder
        elif target[0] == "job":
            folder = target[1].handle.folder if target[1].handle is not None else None
        else:
            folder = target[1].folder
        if folder is None:
            return
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        try:
            subprocess.Popen([opener, str(folder)], stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
        except OSError as exc:
            notify(f"Could not open {folder}: {exc}", level="warn")

    def _reopen(workflow_id):
        """The handle for an offered, referenced or live workflow; ``None`` if it is
        unknown."""
        entry = ctx.recovery.get(workflow_id)
        if entry is not None:
            return entry.handle
        ref = ctx.workflow_refs.get(workflow_id)
        if ref is not None:
            return submission.open_workflow(ref.folder)
        summary = _live_summary(workflow_id)
        return submission.open_workflow(summary["folder"]) if summary is not None else None

    def attach_workflow(workflow_id):
        """Follow a live workflow from here: a new job under the current generation whose
        row can cancel it through the handle and whose completion is delivered once. Live
        means queued or running in the daemon, or a Slurm job that still exists; a record
        without one cannot be attached and nothing is resumed."""
        handle = _reopen(workflow_id)
        if handle is None:
            return None
        job = jobs.find_workflow(workflow_id)
        if job is not None:
            return job
        try:
            state = handle.status
        except OSError as exc:
            notify(f"Workflow {workflow_id} cannot be read: {exc}", level="danger")
            return None
        if handle.done or state not in (QUEUED, RUNNING):
            notify(f"{handle.label} is {state}; it cannot be attached", level="warn")
            return None
        try:
            deliver = _import_result(handle)
        except ValueError as exc:
            notify(f"{handle.label} cannot be attached: {exc}", level="warn")
            return None
        job = _admit(handle.label, heavy=True, detached=True)
        job.handle = handle
        job.started = handle.to_dict().get("started")
        job.token.add_callback(handle.cancel)
        ctx.workflow_refs.attach(handle, handle.label)
        ctx.recovery.forget(workflow_id)
        pulse(job_modal, "close")
        _subscribe(job, handle, deliver)
        _follow_progress(job, handle, _import_progress(handle))
        update_all(canvas)
        return job

    def import_workflow(workflow_id):
        """Load a finished workflow's saved result off the document thread and add its
        spectra once; a workflow already imported is refused."""
        handle = _reopen(workflow_id)
        if handle is None:
            return None
        if ctx.workflow_refs.imported(workflow_id):
            notify(f"{handle.label} was already imported", level="warn")
            ctx.recovery.forget(workflow_id)
            update_all(canvas)
            return None
        try:
            deliver = _import_result(handle)
        except ValueError as exc:
            notify(f"{handle.label} cannot be imported: {exc}", level="warn")
            return None

        def done(result):
            deliver(result)
            ctx.workflow_refs.attach(handle, handle.label, imported=True)

        ctx.recovery.forget(workflow_id)
        recovery_bridge.value = processes.recovery_html(ctx.recovery)
        return load_async(f"import {handle.label}", lambda token=None: handle.result(), done)

    def forget_workflow(workflow_id):
        if ctx.recovery.forget(workflow_id) is not None:
            recovery_bridge.value = processes.recovery_html(ctx.recovery)
            update_all(canvas)

    def _reattach(refs):
        """Adopt a loaded slot's workflow references: a live one is followed again (a job
        already following it is re-subscribed under the new generation), a finished one
        not yet imported is offered for import, an interrupted, failed or cancelled one is
        shown as such. Returns the outcome lines for the load report and the references
        whose record is gone."""
        ctx.workflow_refs.replace(refs)
        outcomes, gone = [], []
        for ref in ctx.workflow_refs:
            handle = submission.open_workflow(ref.folder)
            try:
                state = handle.status
            except (OSError, ValueError) as exc:
                gone.append(session.Skipped("workflow", ref.id, ref.name, ref.folder,
                                            "no record" if isinstance(exc, FileNotFoundError)
                                            else str(exc)))
                continue
            job = jobs.find_workflow(ref.id)
            if job is not None and not handle.done:
                job.generation = ctx.generation      # the slot names it: subscribe anew
                ctx.recovery.forget(ref.id)
                outcomes.append(f"{ref.name} ({state}, followed)")
            elif not handle.done and state in (QUEUED, RUNNING):
                attach_workflow(ref.id)
                outcomes.append(f"{ref.name} ({state}, reattached)")
            elif ref.imported:
                continue                              # its spectra are in the slot itself
            elif _offer(handle, referenced=True) is not None:
                outcomes.append(f"{ref.name} ({state}"
                                + (", result to import)" if state == FINISHED else ")"))
        return outcomes, gone

    ctx.discover_workflows = discover_workflows
    ctx.attach_workflow = attach_workflow
    ctx.import_workflow = import_workflow

    def session_closed(session_context):
        autosave(now=True)      # the document's last state, for Resume last session
        ctx.closed = True
        # Closing a tab cancels this document's work (loads, fits, analysis jobs), except
        # the jobs that asked not to be: every synthesis is a durable workflow that keeps
        # running and writes its record and result under its `_workflows/` folder, where
        # a later document can reopen it by id.
        for job in jobs.active:
            if not job.detached:
                jobs.request_cancel(job.pid)
        if live_view["unwatch"] is not None:
            live_view["unwatch"]()

    ctx.doc.on_session_destroyed(session_closed)

    def _settle(job, result, on_done, exc, on_error=None):
        """Deliver a job's outcome on the document thread. Called from the worker (or,
        for a queued cancellation, the cancelling thread); the state and any delivery
        are applied on the document, guarded by its lifetime and the workspace
        generation the job was admitted under. A closed document only records."""
        def apply():
            if job.token.cancelled:
                jobs.settle(job.pid, CANCELLED)   # user cancelled: drop whatever arrived
            elif exc is not None:
                logs = job.handle.logs if job.handle is not None else ()
                jobs.settle(job.pid, FAILED, error=f"{type(exc).__name__}: {exc}", logs=logs)
            else:
                jobs.settle(job.pid, FINISHED)
            if ctx.closed:
                return                       # recorded; there is no document to tell
            if job.state == FAILED:
                notify(f"{job.name} failed: {exc} — see Processes", level="danger")
                if on_error is not None and job.generation is ctx.generation:
                    try:
                        on_error(exc)
                    except Exception:
                        log.exception("%s (job %d): failure not applied", job.name, job.pid)
            elif job.state == FINISHED:
                if job.generation is not ctx.generation:
                    log.info("%s (job %d) finished after its workspace was replaced; dropped",
                             job.name, job.pid)
                elif on_done is not None:
                    try:
                        on_done(result)      # add_object -> on_change -> update_all
                    except Exception as err:
                        log.exception("Failed to add %s (job %d)", job.name, job.pid)
                        notify(f"{job.name} could not be added: {err}", level="danger")
            update_all(canvas)               # clear the spinner (and show the new object)

        if ctx.closed:
            apply()                      # no document to deliver to: record only
            return
        try:
            ctx.doc.add_next_tick_callback(apply)
        except Exception:
            apply()                      # the document went away under us

    def cancel_pending(pid):
        """Cancel a job from its row in the Processes panel. A still-queued job settles
        at once; a running one is told through its token (a dispatch child is killed, a
        shared dispatch job continues for its other subscribers) and its row shows
        the request until the worker settles. Its result, if any still arrives, is
        dropped in _settle."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        job = jobs.request_cancel(pid)
        if job is None:
            return
        update_all(canvas)
        if job.terminal:
            notify(f"Cancelled {job.name}", level="info")
        else:
            notify(f"Cancelling {job.name}…", level="info")

    def dismiss_job(pid):
        """Drop a retained (failed) row from the Processes panel."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        if jobs.dismiss(pid) is not None:
            update_all(canvas)

    def on_select(filepath):
        remember_recent(filepath)
        name = Path(filepath).name
        def focus_state():
            return (ctx.focus_kind, canvas.selected_uid, tools.active_name(),
                    canvas.main.x_range.start, canvas.main.x_range.end)
        initial = focus_state()
        def add_loaded(spec):
            select = focus_state() == initial
            canvas.add_object(spec, name=name, source=("file", str(filepath)),
                              meta=getattr(spec, "meta", None), select=select)
            if not select:
                notify(f"Opened {name}; select it in Objects to view it.", level="info")
        load_async(name, lambda token=None: load_file(filepath), add_loaded, path=filepath)

    # client -> server channel for the sidebar hide/remove buttons ("action:uid:nonce")
    actions = B.action("object-actions")

    def set_color(value):
        obj = canvas.selected()
        if obj is not None:
            canvas.set_style(obj, color=value)

    # clicking a line-position indicator opens the line list in the bottom half of the
    # workspace, scrolled to that transition's row (lineannot -> "line-picked" -> linenav)
    ctx.subscribe("line-picked",
                  lambda lam0_aa, ion=None, **_: linenav.open_at(lam0_aa, ion))

    # action name -> handler(value). The mutually-exclusive sidebar tools go through the
    # ToolManager, which owns "activate one, stand down the rest" (see toolmanager.py).
    ACTIONS = {
        "toggle":       canvas.toggle_object,
        "remove":       canvas.remove_object,
        "cancel-pending": cancel_pending,
        "dismiss-job":  dismiss_job,
        "workflow-discover": lambda v: discover_workflows(),
        "workflow-attach": attach_workflow,
        "workflow-import": import_workflow,
        "workflow-forget": forget_workflow,
        "workflow-cancel": cancel_workflow,
        "job-details":  show_details,
        "folder-open":  open_folder,
        "run-cancel":   cancel_run,
        "run-departures": use_run_departures,
        "select":       focus_spectrum,
        "color":        set_color,
        "profile":      lambda v: range_tool.ew_setting(range_tool.ew.set_profile, v),
        "ew-snr":       lambda v: range_tool.ew_setting(range_tool.ew.set_snr, v),
        "ew-cont":      lambda v: range_tool.ew_setting(range_tool.ew.set_cont, v),
        "range-ew":     lambda v: range_tool.toggle_ew(),
        "mask-toggle":  lambda v: tools.toggle("range"),
        "range-cancel": lambda v: range_tool.cancel(),
        "range-new-mask": lambda v: range_tool.new_mask(),
        "range-synth":  lambda v: range_tool.synthesise(),
        "range-fit":    lambda v: range_tool.fit(),
        "range-add":    lambda v: range_tool.add_window(),
        "synth-new":    lambda v: synth_modal.open_range(),
        "synth-submit": lambda v: synth_modal.submit(),
        "annot-toggle": lambda v: ctx.annot.toggle_visible(),
        "linejump":     linenav.jump,
        "linesort":     linenav.set_sort,
        "lineclose":    lambda v: linenav.close(),
        "lineopen":     lambda v: linenav.close() if linenav.active else linenav.open_at(),
        "loadrun":      lambda v: load_current_run(),
        "filter":       lambda v: tools.toggle("filter"),
        "filter-clear": lambda v: filtertool.clear(),
        "filtertoggle": lambda v: filtertool.set_visible(v == "on"),
        "normalise":    lambda v: tools.toggle("normalise"),
        "norm-init":    lambda v: normtool.initialise(),
        "norm-clear":   lambda v: normtool.clear(),
        "linemask":     lambda v: tools.toggle("linemask"),
        "abundance":    lambda v: tools.toggle("abundance"),
        "abundcorr":    lambda v: tools.toggle("abundcorr"),
        "atom":         lambda v: tools.toggle("atom"),
        "ionbal":       lambda v: tools.toggle("ionbal"),
    }

    def on_action(attr, old, new):
        if not new:
            return
        action, value = new.split(":")[:2]
        handler = ACTIONS.get(action)
        if handler is not None:
            handler(value)

    actions.on_change("value", on_action)

    # --- session save slots (named saves kept in the config file) ---------
    # a dedicated space-delimited channel: names may contain spaces/colons, so the
    # ":"-delimited object-actions bridge can't carry them safely (cf. cache-actions).
    session_actions = B.action("session-actions")
    # broadcast to every session list on the page: the Sessions modal *and* the inline
    # empty-state loader shown over the canvas when nothing is loaded
    sessions_bridge = B.html_all(".js-session-results", name="session-render")
    sessions_bridge.value = sessions_html(store.list_sessions())  # seed initial list

    def render_sessions():
        """Re-read the slots into every session list; a config the store rejects leaves
        the list as it was and reports why."""
        try:
            sessions_bridge.value = sessions_html(store.list_sessions())
        except store.ConfigError as exc:
            logging.exception("Could not list the saved sessions")
            notify(f"Could not read saved sessions: {exc}", level="danger")

    saved_workspace = {"snapshot": snapshot(), "pending": None}
    load_prompt = B.dialog("workspace-load-modal")
    sessions_dialog = B.dialog("session-modal")

    # --- autosave: the live workspace kept in one slot, so a closed tab resumes ---------
    # Every update_all (a load, a result, a view change) schedules a write once the burst
    # settles; the session-destroyed hook writes at once. A blank workspace is never
    # written, so a fresh tab does not replace the last real session, and a snapshot equal
    # to the last one written is skipped. Failures are logged: this must never interrupt
    # the operation that triggered it.
    AUTOSAVE_DELAY_MS = 2000
    autosave_delay = Trailing(AUTOSAVE_DELAY_MS)
    autosaved = {"snapshot": None}

    def write_autosave():
        try:
            snap = snapshot()
            if session.is_blank(snap) or snap == autosaved["snapshot"]:
                return
            store.save_autosave(snap, workflows=ctx.workflow_refs.to_list())
            autosaved["snapshot"] = snap
        except Exception:
            logging.exception("Autosave failed")

    def autosave(now=False):
        if now:
            autosave_delay.cancel("autosave")
            write_autosave()
        else:
            autosave_delay("autosave", write_autosave)

    def on_session_action(attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]  # strip the trailing nonce
        # every store call may raise ConfigError (unreadable/invalid config, failed
        # replacement, lock timeout); the message names the path and the cause, and the
        # live workspace stays as it is so the user can retry after fixing the cause.
        if kind == "save":
            name = arg.strip() or time.strftime("Save %Y-%m-%d %H:%M")
            try:
                store.save_session(name, snapshot(), workflows=ctx.workflow_refs.to_list())
            except store.ConfigError as exc:
                logging.exception("Could not save the session %r", name)
                notify(f"Save “{name}” failed: {exc}", level="danger")
                return
            saved_workspace["snapshot"] = snapshot()
            render_sessions()
            notify(f"Saved “{name}”", level="success")
        elif kind in ("load", "load-discard", "load-save"):
            if kind == "load" and snapshot() != saved_workspace["snapshot"]:
                saved_workspace["pending"] = arg
                pulse(load_prompt, "open")
                return
            if kind != "load":
                arg = saved_workspace["pending"]
                if arg is None:
                    return
                if kind == "load-save":
                    try:
                        store.save_session(time.strftime("Workspace %Y-%m-%d %H:%M"),
                                           snapshot(), workflows=ctx.workflow_refs.to_list())
                    except store.ConfigError as exc:
                        notify(f"Could not save workspace: {exc}", level="danger")
                        return
                pulse(load_prompt, "close")
            try:
                slot = store.get_slot(int(arg))
            except store.ConfigError as exc:
                logging.exception("Could not read the saved session %s", arg)
                notify(f"Could not read the save: {exc}", level="danger")
                return
            if slot is None:
                notify("That save is no longer available.", level="danger")
                render_sessions()
            else:
                report = restore_workspace(slot["snapshot"], slot["name"], slot["workflows"])
                if report is not None:
                    saved_workspace["snapshot"] = snapshot()
                    saved_workspace["pending"] = None
                    pulse(sessions_dialog, "close")
        elif kind == "delete":
            try:
                found = store.delete_session(int(arg))
            except store.ConfigError as exc:
                logging.exception("Could not delete the saved session %s", arg)
                notify(f"Delete failed: {exc}", level="danger")
                return
            render_sessions()
            if not found:
                notify("That save was already gone.", level="warn")
        elif kind == "refresh":
            render_sessions()  # re-read slots (another tab may have changed them)

    session_actions.on_change("value", on_session_action)

    # separate channel for the cache page, format "kind <arg> <nonce>" (space-delimited):
    # kind/nonce have no spaces, so partition+rsplit recover an arg path even with spaces.
    # The ":"-delimited object-actions bridge can't carry a filesystem path safely.
    cache_actions = B.action("cache-actions")

    def on_cache_action(attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]  # strip the trailing nonce
        if kind == "del":
            synth_modal.delete_cache(arg)
        elif kind == "clear":
            synth_modal.clear_cache()
        elif kind == "refresh":
            synth_modal.refresh_cache()

    cache_actions.on_change("value", on_cache_action)

    # the "Load this run" button appears when the picker is inside a tumult run folder
    def on_picker_dir(folder):
        run_name_bridge.value = html.escape(Path(folder).name) if folder and is_run(folder) else ""
        toggle(run_ready_bridge, bool(folder and is_run(folder)))

    def load_current_run():
        folder = getattr(picker, "current_dir", None)
        if folder is None or not is_run(folder):
            return
        remember_recent(folder)
        name = Path(folder).name

        def add(spec):
            color, dash = canvas.next_run_style()  # picked at add time -> no colour clashes
            canvas.add_object(spec, name=name, color=color, dash=dash,
                                   source=("run", str(folder)))

        load_async(name, lambda token=None: load_run(folder), add)

    # `busy_provider`: on_select hands the read to a worker thread and returns at once, so
    # the row's spinner is carried by the job table until the object lands
    picker = FilePicker(on_select, container_id="file-results", start_dir=start_dir,
                        on_dir=on_picker_dir, recent_provider=store.recent_files,
                        busy_provider=jobs.paths)
    picker.layout.name = "file_picker"

    # the model-atom picker of the Open Model Atom modal, rooted at input/atoms
    atom_picker = FilePicker(on_atom_select, container_id="atom-file-results",
                             start_dir=_atoms_dir if _atoms_dir.is_dir() else Path.home(),
                             busy_provider=jobs.paths)
    atom_picker.layout.name = "atom_file_picker"

    def sync_linelist():
        prefs = load_prefs()
        path = str(_resolve_path(prefs, prefs["linelist"]))
        strengths.set_master_path(path)
        abundancetool.master_path = path
        ionbaltool.master_path = path
        auto = masktool.source("auto")
        if auto is not None:
            auto.linelist_path = path
        ctx.annot.resync()
        linenav.render()

    synth_modal.on_inputs_change = sync_linelist
    if store.get_pref("master_linelist") is not None:
        with store.transaction() as data:
            data["preferences"].pop("master_linelist", None)
        notify("Navigation and analysis now use the synthesis Line list in Preferences.")

    from draft import CalculationDraft
    synth_modal.draft = CalculationDraft(setup=synth_modal.setup, atoms=synth_modal.atoms,
        owner=synth_modal.fit, attributes=("vals", "flags", "abund_rows", "target"))
    abundancetool.draft = CalculationDraft(forms=[("abundance", abundancetool.form1)],
        setup=abundancetool.setup, atoms=abund_atoms, owner=abundancetool, attributes=("ref_uid",))
    abundcorrtool.draft = CalculationDraft(forms=[("abundcorr", abundcorrtool.form1)],
        setup=abundcorrtool.setup)
    ionbaltool.draft = CalculationDraft(forms=[("ionbal", form) for form in
        (ionbaltool.form1, ionbaltool.form2, ionbaltool.formadv, ionbaltool.setup.form)])
    drafts = {"synth-modal": synth_modal.draft, "abund-modal": abundancetool.draft,
              "corr-modal": abundcorrtool.draft, "ionbal-modal": ionbaltool.draft}
    dialog_actions = B.action("dialog-actions")
    def close_draft(attr, old, new):
        draft = drafts.get(new.split(" ")[0])
        if draft is not None:
            draft.discard()
    dialog_actions.on_change("value", close_draft)

    # every html control is a hidden model collected under one root: the bridges/actions
    # `B` gathered as they were created, plus each component's own form/bridge models
    bridges = column(*B.models,
                     *settings.models, *inspector.models, *linenav.models, *abund_atoms.models,
                     *linejump_form.models, *refstars_page.models,
                     *synth_modal.models, *range_tool.models, *normtool.models,
                     *filtertool.models, *masktool.models,
                     *abundancetool.models, *abundcorrtool.models,
                     *ionbaltool.models)
    bridges.name = "bridges"

    doc = curdoc()
    doc.title = "TUMULT Studio"
    for root in (canvas.nav, canvas.main, filtertool.fig,
                 *abundancetool.figs, *abundcorrtool.figs, abundcorrtool.cog_fig,
                 *ionbaltool.figs, picker.layout,
                 *masktool.roots,
                 synth_modal.atoms.picker.layout, abund_atoms.picker.layout,
                 *(p.layout for p in synth_modal.pref_pickers.values()),
                 synth_modal.grids.picker.layout, synth_modal.atmos_picker.layout,
                 abundancetool.setup.atmos_picker.layout,
                 abundcorrtool.setup.atmos_picker.layout,
                 atom_picker.layout, bridges):
        doc.add_root(root)

    ionbaltool.refresh()      # the balance modal's target selector

    repaint_processes()       # the daemon's jobs are rows from the first paint on

    inspector.bind(canvas.selected())
    sel = canvas.selected()
    active_color = sel.line.glyph.line_color if sel else None
    _auto = masktool.source("auto")     # None while masktool.AUTO_MASKS is off
    doc.template_variables.update(
        css=(HERE / "static" / "style.css").read_text(),
        main_empty=not canvas.objects,  # bake the empty-state class into the first paint
        objects_html=objects_panel_html(canvas),
        processes_html=processes_bridge.value,
        processes_usage_html=processes_usage_bridge.value,
        inspector_html=inspector_html(canvas),
        footer_html=footer_html(canvas),
        swatches_html=swatch_html(active_color),
        edit_html=inspector.edit_html(),
        style_html=inspector.style_form.html(),
        settings_html=settings.html(),
        synth_run_html=synth_modal.run_form.html(),
        synth_context_html=synth_modal.context_html(),
        synth_actions_html=synth_modal.actions_html(),
        synth_abund_html=synth_modal.abund_form.html(),
        synth_mode_html=synth_modal.mode_form.html(),
        synth_prefs_html=synth_modal.pref_forms["synthprefs"].html(),
        synth_engine_html=synth_modal.pref_forms["synthengine"].html(),
        pref_picker_results=synth_modal.pref_picker_results(),
        pref_labels=PREF_LABELS,
        pref_missing=set(synth_modal.missing),
        pref_show_slurm=load_prefs().get("workflow_runner") == "slurm",
        grid_list_html=synth_modal.grids.list_html(),
        refstars_html=refstars_page.list_html(),
        grid_results_name=synth_modal.grids.results_name,
        synth_grid_html=synth_modal.grid_ui_html(),
        synth_show_file=synth_modal.setup.show_file(),
        synth_file_results_name=synth_modal.atmos_picker.results_name,
        synth_cache_html=synth_modal.cache_html(),
        synth_warning=synth_modal.has_warning,
        synth_atoms_html=synth_modal.atoms.chips_html(),
        synth_is_nlte=synth_modal.is_nlte,
        atom_results_name=synth_modal.atoms.results_name,
        fit_html=synth_modal.fit.html(),
        annot_label=ctx.annot.button_html(),
        linesettings_html=linejump_form.html(),
        linemask_results_name=masktool.source("file").results_name,
        # the Auto-generate pane is only rendered while that source exists
        # (masktool.AUTO_MASKS); its fields are undefined otherwise
        linemask_auto=_auto is not None,
        linemask_auto_results_name=_auto.results_name if _auto else "",
        linemask_auto_ref_html=_auto.ref_html() if _auto else "",
        linemask_auto_path=html.escape(Path(_auto.linelist_path).name) if _auto else "",
        linemask_auto_html=_auto.form.html() if _auto else "",
        linemask_objects_html=masktool.objects_html(),
        sessions_html=sessions_html(store.list_sessions()),
        abund_results_html=abundancetool.results_html(),
        abund_table_html="",
        abund_ref_html=abundancetool.ref_html(),
        abund1_html=abundancetool.form1.html(),
        abund_context_html="",
        abund_actions_html=abundancetool.actions_html(),
        abund_setup_note_html=abundancetool.setup_note_html(),
        abund_params_html=abund_setup.run_form.html(),
        abund_base_html=abund_setup.base_form.html(),
        abund_mode_html=abund_setup.mode_form.html(),
        abund_grid_html=abund_setup.grid_ui_html(),
        abund_show_file=abund_setup.show_file(),
        abund_atoms_html=abund_atoms.chips_html(),
        abund_atoms_results_name=abund_atoms.results_name,
        abund_is_nlte=abund_setup.is_nlte,
        abund_file_results_name=abund_setup.atmos_picker.results_name,
        atom_objects_html=atomtool.objects_html(),
        atom_file_results_name=atom_picker.results_name,
        corr_table_html="",
        corr1_html=abundcorrtool.form1.html(),
        corr_params_html=corr_setup.run_form.html(),
        corr_grid_html=corr_setup.grid_ui_html(),
        corr_show_file=corr_setup.show_file(),
        corr_file_results_name=corr_setup.atmos_picker.results_name,
        ionbal_results_html=ionbaltool.results_html(),
        ionbal_table_html="",
        ionbal_context_html="",
        ionbal1_html=ionbaltool.form1.html(),
        ionbal2_html=ionbaltool.form2.html(),
        ionbalsetup_html=ionbaltool.setup.form.html(),
        ionbaladv_html=ionbaltool.formadv.html(),
        form_script=FORM_SCRIPT + NUM_INPUT_SCRIPT,
        atmos_script=ATMOS_SCRIPT,
        analysis_script=ANALYSIS_SCRIPT,
        results_name=picker.results_name,
    )
    return ctx        # local harnesses drive the document through its context
