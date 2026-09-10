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

import numpy as np  # noqa: E402
from bokeh.io import curdoc  # noqa: E402
from bokeh.layouts import column  # noqa: E402

from tumult.spectrum_io import Spectrum, load_spectrum, load_file, is_run, load_run  # noqa: E402
from canvas import SpectrumCanvas, DASH_CSS, LINE  # noqa: E402
from context import AppContext  # noqa: E402
from controls import Field, Form, FORM_SCRIPT, NUM_INPUT_SCRIPT  # noqa: E402
from inspector import Inspector, swatch_html  # noqa: E402
from ew import EquivalentWidth  # noqa: E402
from filtertool import FilterTool  # noqa: E402
from normalisetool import NormaliseTool  # noqa: E402
from maskstore import MaskStore  # noqa: E402
from masktool import MaskTool  # noqa: E402
from linenav import LineNavigator  # noqa: E402
from lineannot import LineAnnotations  # noqa: E402
from tumult import refstars  # noqa: E402
from tumult.linestrength import LineStrengths  # noqa: E402
from synthpanel import SynthPanel  # noqa: E402
from synthsetup import SynthSetup  # noqa: E402
from abundancetool import AbundanceTool  # noqa: E402
from abundwizard import PROFILES  # noqa: E402  (shared EW-profile choices)
from abundcorrtool import AbundCorrTool  # noqa: E402
from ionbalance import IonBalanceTool  # noqa: E402
from tumult.synthesis import workflows as synth  # noqa: E402  (CancelToken/CancelledError for cancellable loads)
from filepicker import FilePicker  # noqa: E402
from toolmanager import ToolManager  # noqa: E402
from session import Workspace, ObjectSnap, FilterSnap  # noqa: E402
from htmlbridge import BridgeSet, pulse, toggle  # noqa: E402
from debounce import Trailing  # noqa: E402
from units import X_UNITS, Y_UNITS, to_range  # noqa: E402
from tumult import linestrength
from tumult.synthesis.runtime import IO as _EXEC, WORKFLOWS as _SYNTH_EXEC, HOP as _HOP_EXEC
from profiling import profiled  # noqa: E402
import store  # noqa: E402

# how long the view has to hold still before the visible band's line strengths are computed
STRENGTH_DELAY_MS = 400

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(name)s: %(message)s")


class FooterErrorHandler(logging.Handler):
    """Surface any logged error (a caught exception, an ``logging.exception``/``error``
    call) as a footer warning that points the user at the console for the traceback.

    The console still gets the full record via the root logger's stream handler; this
    just flags that *something* went wrong so it isn't missed. The notify is scheduled
    onto the document thread, so it is safe to fire from background worker threads."""

    def __init__(self, doc, notify):
        super().__init__(level=logging.ERROR)
        self._doc = doc
        self._notify = notify

    def emit(self, record):
        def show():
            try:
                self._notify("Something went wrong — check the console for details.",
                             level="danger")
            except Exception:
                pass  # never let the footer notify re-raise into the logging machinery
        try:
            self._doc.add_next_tick_callback(show)
        except Exception:
            pass  # document/session already gone, or its callback queue is closed


def picker_start_dir():
    """Where the file picker opens (the canvas itself starts empty): the folder of
    $TUMULT_SPECTRUM if set, else the most recently opened file's folder, else home."""
    env = os.environ.get("TUMULT_SPECTRUM")
    for candidate in ([env] if env else []) + store.recent_files():
        folder = Path(candidate).parent
        if folder.is_dir():
            return folder
    return Path.home()


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
        n = len((s.get("snapshot") or {}).get("objects") or [])
        try:
            when = datetime.fromtimestamp(s.get("saved_at", 0)).strftime("%Y-%m-%d %H:%M")
        except (OSError, ValueError, OverflowError):
            when = ""
        sub = f'{n} object{"" if n == 1 else "s"}' + (f' · {when}' if when else "")
        rows.append(
            f'<div class="session-item" data-sid="{s["id"]}" title="Load this save">'
            f'<i class="ti ti-device-floppy"></i>'
            f'<div class="session-meta"><span class="session-name">{html.escape(s["name"])}</span>'
            f'<span class="session-sub">{html.escape(sub)}</span></div>'
            f'<button class="session-del" data-sid="{s["id"]}" '
            f'title="Delete this save"><i class="ti ti-trash"></i></button>'
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
            f'title="Remove"><i class="ti ti-x"></i></button>'
            f'</div>'
        )
    return "".join(rows)


# shown in the Objects panel while it holds nothing (and nothing is loading): how to
# get a first spectrum onto the canvas
OBJECTS_EMPTY_HINT = (
    '<div class="object-empty">No spectra yet. Open one with '
    '<b>Open File</b> <i class="ti ti-folder-open"></i> above, or click '
    '<b>Synth</b> to compute a new spectrum.</div>'
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
    return (
        _field("Name", obj.name)
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


def processes_html(pending, heavy=()):
    """The Processes panel body: one row per background run, each with a cancel button
    (data-uid = the pending id) that drops the run and discards its result.

    A run in ``heavy`` is one whose loss actually costs something (a synthesis, a fit),
    and its row asks for confirmation before cancelling; the cheap ones (a file load, a
    line-strength band) cancel on the first click."""
    return "".join(
        f'<div class="object-item process" data-pid="{pid}">'
        '<span class="spinner"></span>'
        f'<span class="object-name" title="{html.escape(name)}">{html.escape(name)}</span>'
        f'<button class="object-btn danger" data-action="cancel-pending" data-uid="{pid}" '
        f'data-confirm="{1 if pid in heavy else 0}" '
        f'title="Cancel"><i class="ti ti-x"></i></button>'
        f'</div>'
        for pid, name in pending.items()
    )


def build():
    from linenav import MASTER_PATH

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
    processes_bridge = B.html("processes")            # the running-work list, bottom left
    processes_show_bridge = B.cls(".processes-panel", "show")  # hidden while nothing runs
    inspector_bridge = B.html("inspector")
    footer_bridge = B.html("footer-info")
    status_bridge = B.html("footer-status")  # transient info messages (left of footer)

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
            curdoc().add_timeout_callback(restore, revert_ms)
        except Exception:
            pass

    # route any logged error to the footer so caught exceptions aren't lost to the
    # console alone; tear it down with the session so dead documents aren't notified.
    # Capture the logger in a local (the closure must not rely on the module global
    # `logging`, which Bokeh has already cleared by session-destroy time).
    _root_logger = logging.getLogger()
    _err_handler = FooterErrorHandler(curdoc(), notify)
    _root_logger.addHandler(_err_handler)
    curdoc().on_session_destroyed(lambda ctx: _root_logger.removeHandler(_err_handler))

    # (the EW / Filter / Normalise tools create and own their bridges via the context)
    linelist_bridge = B.html("linelist-table")
    linelist_mode_bridge = B.cls("#main", "linelist-mode")  # the line-list bottom panel
    run_name_bridge = B.html("run-name")          # name shown on the "Load run" button
    run_ready_bridge = B.cls("#load-run", "ready")  # reveal it inside a run folder
    loading_bridge = B.cls("#footer", "loading")       # footer progress bar while loading
    empty_bridge = B.cls("#main", "empty")  # swap the empty canvas for the inline session loader
    synth_mode_bridge = B.cls(".sidebar-right", "synth-active")  # swap to the synth panel
    synth_btn_bridge = B.cls("#synth-btn", "active")             # toolbar button state
    synth_warn_bridge = B.cls("#settings-warning", "show")       # missing-files sign on the cog
    # (the abundance / correction wizards create and own their bridges via the context)

    # which left-sidebar list owns the single selection ("spectrum" | "result" | "corr");
    # the wizards drive it through ctx.focus_kind + the "focus-changed" event, and the
    # Objects list repaints to follow the highlight.
    ctx.focus_kind = "spectrum"
    ctx.subscribe("focus-changed",
                  lambda kind, **_: setattr(objects_bridge, "value",
                                            objects_panel_html(canvas, kind)))

    def synth_on_mode(state):
        pulse(synth_mode_bridge, state)
        pulse(synth_btn_bridge, state)

    def synth_on_warn(has_warning):
        toggle(synth_warn_bridge, has_warning)

    inspector = linenav = abundancetool = abundcorrtool = ionbaltool = None  # assigned after the canvas
    warm_strengths = None   # ditto: update_all may run before it is defined
    picker = None           # ditto: the spectrum file picker, whose rows spin while loading
    pending = ctx.pending  # id -> name of each background run (the Processes panel)
    pending_futures = {}  # id -> Future, so a still-queued run can be cancelled before it starts
    pending_tokens = {}   # id -> synth.CancelToken, to hard-kill a running synthesis
    # id -> the file a pending load came from, for the picker's per-row spinners. Separate
    # from `pending` because that one is a display name and a synthesis has no file at all.
    pending_paths = {}
    pending_heavy = set()  # ids whose cancel button confirms first (synthesis, fits)

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
        return [r.start, r.end] if r.start is not None and r.end is not None else None

    def snapshot():
        objs = []
        for o in canvas.objects:
            if not o.source:
                continue
            kind = o.source[0]
            # synth runs are interpolated (no backing file), so embed their data;
            # file/run objects need a path to reload from. Anything else can't be restored.
            if kind == "synth":
                data = o.raw.to_dict()
            elif o.source[1]:
                data = None
            else:
                continue
            g = o.line.glyph
            knots = ({"x": list(map(float, o.norm_knots[0])),
                      "y": list(map(float, o.norm_knots[1]))}
                     if o.norm_knots is not None else None)
            # synth params are tiny and have no source file, so persist them; FITS metadata
            # is large and re-read from disk on restore, so isn't stored.
            meta = o.meta if kind == "synth" else None
            objs.append(ObjectSnap(uid=o.uid,
                kind=kind, path=o.source[1] or "", name=o.name, data=data, meta=meta,
                color=g.line_color, dash=o.dash, width=g.line_width, alpha=g.line_alpha,
                visible=o.visible, conv=o.conv, norm_knots=knots,
                already_normalised=o.already_normalised))
        sel = canvas.selected()
        idx = next((i for i, o in enumerate(canvas.objects) if o is sel), 0)
        configured = filtertool is not None and filtertool._configured
        # only persist a low-pass cutoff when a low-pass filter actually exists
        has_lowpass = (filtertool is not None and filtertool.method == "lowpass"
                       and any(o.filtered is not None for o in canvas.objects))
        return Workspace(
            objects=objs,
            units={"x": canvas.xlabel, "y": canvas.ylabel},
            view={"main": _range(canvas.main.x_range), "nav": _range(canvas.nav.x_range)},
            selected=idx,
            filter=FilterSnap(
                active=canvas.show_filtered,
                method=filtertool.method if filtertool else "brault",
                res=filtertool.form.value("res") if has_lowpass else None,
                fmax=filtertool.fmax if configured else None,
                floor=filtertool.floor if configured else None),
            abundance=abundancetool.snapshot() if abundancetool is not None else {},
            abundcorr=abundcorrtool.snapshot() if abundcorrtool is not None else {},
            ionbal=ionbaltool.snapshot() if ionbaltool is not None else {},
        ).to_dict()

    def update_all(canvas):
        objects_bridge.value = objects_panel_html(canvas, ctx.focus_kind)
        processes_bridge.value = processes_html(pending, pending_heavy)
        toggle(processes_show_bridge, bool(pending))
        inspector_bridge.value = inspector_html(canvas)
        footer_bridge.value = footer_html(canvas)
        synth_panel.refresh_targets()  # keep the Fit-tab target dropdown in sync
        toggle(loading_bridge, bool(pending))
        # show the inline session loader while the canvas holds nothing (and nothing
        # is on its way in)
        toggle(empty_bridge, not canvas.objects and not pending)
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

    inspector = Inspector(canvas, on_fit=edit_fit, notify=notify)
    ew = EquivalentWidth(ctx)

    master = store.get_pref("master_linelist")
    if not master or not Path(master).is_file():
        master = MASTER_PATH
    # where every strength comes from: HOP, on the MARCS model nearest the object's star,
    # computed per band of spectrum on demand (see linestrength for why it matters)
    strengths = LineStrengths(
        master_path=master,
        refstar=store.get_pref("strength_refstar", refstars.AUTO_KEY))
    ctx.strengths = strengths

    linenav = LineNavigator(
        canvas, strengths, on_table=lambda h: setattr(linelist_bridge, "value", h),
        on_open=lambda on: toggle(linelist_mode_bridge, on))

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
        obj = canvas.selected()
        lo, hi = ctx.annot._window_aa()
        if not strengths.pending(obj, lo, hi):
            return False

        # Snapshot the reference and band on the document thread. Do not read
        # canvas objects or form values from the HOP worker.
        star = strengths.star_for(obj)
        linelist = strengths.linelist
        band = linestrength.band_bounds(lo, hi)
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
        # the mask view's analysis actions: hand the open set to the abundance / balance
        # tools, which are built further down (the lambdas resolve them at click time)
        on_fit_balance=lambda mset: ionbaltool.open_for(mset),
        on_measure_abund=lambda mset: abundancetool.open_for(mset))
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
        synth_panel.refresh_cache()  # the new run is now in the cache listing

    synth_panel = SynthPanel(
        canvas,
        submit=lambda name, loader, on_done: load_async(name, loader, on_done,
                                                        executor=_SYNTH_EXEC),
        add=add_synth,
        on_mode=synth_on_mode,
        on_warn=synth_on_warn,
        notify=notify)

    # the Abundance tool owns its synthesis setup (its own Setup sub-tab), sharing only the
    # app-global grid registry and the Synthesize tab's NLTE model-atom selection; it inherits
    # the Synthesize-tab configuration on first open and is independent thereafter.
    # its fields live in the abundance modal, so every edit restates the setup in the object
    # panel, where the Run synthesis button stands
    _abund_note = lambda *_: (abundancetool is not None and abundancetool.push_setup_note())
    abund_setup = SynthSetup(
        "abund", "abund_run", synth_panel.grids,
        atoms_provider=lambda: list(synth_panel.atoms.selected), notify=notify,
        on_run_change=_abund_note, on_base_change=_abund_note,
        on_mode_change=_abund_note, on_atmos_change=_abund_note)
    abundancetool = AbundanceTool(
        ctx,
        setup=abund_setup,
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        # follow the canvas's Filtered switch: when filtering is on, measure/fit the filtered
        # spectrum (so EW hatches + line fits match what is shown), else the raw observation
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master)

    # the abundance-correction tool: its own synthesis setup (target/NLTE atmosphere) and its
    # own single model-atom selector, sharing only the app-global grid registry (one grid
    # flagged as the 1D LTE reference).
    from atompicker import AtomPicker  # noqa: E402  (second, namespaced instance)
    _corr_atoms_dir = Path(str(store.get_pref("dispatch_cwd") or ".")).expanduser() / "input" / "atoms"
    corr_atoms = AtomPicker(
        start_dir=_corr_atoms_dir if _corr_atoms_dir.is_dir() else Path.home(),
        ns="corr-atoms", section="abundcorr_run", key="atoms",
        picker_sel="#corr-atom-picker", section_sel="#corr-atoms",
        actions_name="corr-atom-actions", results_id="corr-atom-results",
        layout_name="corr_atom_picker", del_class="corr-atom-del", single=True)
    corr_setup = SynthSetup(
        "corr", "abundcorr_run", synth_panel.grids,
        atoms_provider=lambda: list(corr_atoms.selected), notify=notify)

    def add_corr_spectrum(spec, name):
        color, dash = canvas.next_run_style()
        canvas.add_object(spec, name=name, color=color, dash=dash, source=("synth", None))

    abundcorrtool = AbundCorrTool(
        ctx,
        setup=corr_setup,
        grids=synth_panel.grids,
        atom_provider=lambda: list(corr_atoms.selected),
        add_spectrum=add_corr_spectrum)

    # The ionisation-balance modal owns its grid selection and always interpolates it.  Only
    # the app-global registered-grid list and synthesis preferences are shared.
    ionbaltool = IonBalanceTool(
        ctx,
        grids=synth_panel.grids,
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master)

    # single-selection coordination between the Objects and the wizards' result lists
    def focus_spectrum(uid):
        if tools.active_name() in ("abundance", "abundcorr", "ionbal", "linemask"):
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
    tools.register("ew", ew)
    tools.register("synth", synth_panel)
    tools.register("normalise", normtool)
    tools.register("linemask", masktool)
    tools.register("filter", filtertool, runner=run_busy)
    tools.register("abundance", abundancetool)
    tools.register("abundcorr", abundcorrtool)
    tools.register("ionbal", ionbaltool)

    # --- session save / restore -------------------------------------------
    def _restore_range(rng, vals):
        if not vals or vals[0] is None or vals[1] is None:
            return
        lo, hi = vals
        if rng.bounds:
            lo, hi = max(lo, min(rng.bounds)), min(hi, max(rng.bounds))
        if lo < hi:
            rng.start, rng.end = lo, hi

    def restore_workspace(snap):
        ws = Workspace.from_dict(snap)
        if not ws.objects:
            return
        filtertool.clear()
        for o in list(canvas.objects):
            canvas.remove_object(o.uid)
        for obj in ws.objects:
            try:
                spec = _load_obj(obj)
            except Exception:
                logging.exception("Could not restore %s", obj.path)
                continue
            # synth meta is restored from the slot; FITS meta rides along on the reloaded spec
            meta = obj.meta if obj.kind == "synth" else getattr(spec, "meta", None)
            o = canvas.add_object(spec, name=obj.name, color=obj.color, dash=obj.dash,
                                       source=(obj.kind, obj.path or None), meta=meta, uid=obj.uid)
            if obj.width is not None or obj.alpha is not None:
                canvas.set_style(o, width=obj.width, alpha=obj.alpha)
            knots = obj.norm_knots
            o.already_normalised = bool(obj.already_normalised)
            if knots and knots.get("x"):
                o.norm_knots = (np.asarray(knots["x"], float),
                                np.asarray(knots["y"], float))
            if o.already_normalised or (knots and knots.get("x")):
                canvas.normalise_object(o)
            o.conv = obj.conv or {}
            if o.conv or o.already_normalised or (knots and knots.get("x")):
                canvas._refresh_object(o)  # apply conv/normalisation to main + nav
            if not obj.visible:
                canvas.toggle_object(o.uid)
        if not canvas.objects:
            return
        u = ws.units
        if u.get("x"):
            canvas.set_xunit(u["x"])
        if u.get("y"):
            canvas.set_yunit(u["y"])
        _restore_range(canvas.nav.x_range, ws.view.get("nav"))
        _restore_range(canvas.main.x_range, ws.view.get("main"))
        idx = ws.selected
        if 0 <= idx < len(canvas.objects):
            canvas.select_object(canvas.objects[idx].uid)
        f = ws.filter
        if f.method == "lowpass" and f.res is not None:
            filtertool.apply_saved_lowpass(f.res)
            if not f.active:
                filtertool.set_visible(False)
        elif f.fmax is not None and f.floor is not None:
            filtertool.apply_saved(f.fmax, f.floor)
            if not f.active:
                filtertool.set_visible(False)
        settings.set_values(xunit=canvas.xlabel, yunit=canvas.ylabel)  # sync the prefs form
        # restore the wizard results last: the abund/params tools re-link each result's
        # reference to the just-restored observations (by name), so objects must exist first.
        if abundancetool is not None:
            abundancetool.restore(ws.abundance)
        if abundcorrtool is not None:
            abundcorrtool.restore(ws.abundcorr)
        if ionbaltool is not None:
            ionbaltool.restore(ws.ionbal)
        update_all(canvas)

    # --- background loading: placeholder + spinner now, real object when ready ----
    def load_async(name, loader, on_done, executor=_EXEC, path=None):
        pid = _nonce()
        pending[pid] = name
        if executor is _SYNTH_EXEC:
            pending_heavy.add(pid)  # a synthesis or fit: confirm before throwing it away
        if path is not None:
            pending_paths[pid] = str(path)   # a picker row spins for as long as this stands
        token = synth.CancelToken()      # cancel handle: the loader threads it into synth.*
        pending_tokens[pid] = token
        update_all(canvas)               # show the spinner placeholder immediately
        doc = curdoc()

        def work():                      # runs in a worker thread
            if token.cancelled:          # cancelled before the worker even started
                _settle(doc, pid, name, None, None, None)
                return
            try:
                # loaders take the token so synthesis can be hard-cancelled; non-synth
                # loaders accept it and ignore it.
                result = profiled(name, lambda: loader(token))
            except synth.CancelledError:
                _settle(doc, pid, name, None, None, None)  # token.cancelled -> dropped
            except Exception as exc:
                _settle(doc, pid, name, None, None, exc)
            else:
                _settle(doc, pid, name, result, on_done, None)

        pending_futures[pid] = executor.submit(work)
        return pending_futures[pid]

    # Coordinators can overlap; the shared dispatch service enforces the CPU budget.
    ctx.submit = lambda name, loader, on_done: load_async(name, loader, on_done,
                                                          executor=_SYNTH_EXEC)
    ctx.execution = synth.get_runtime()

    def session_closed(session_context):
        ctx.closed = True
        for token in list(pending_tokens.values()):
            token.cancel()
        for future in list(pending_futures.values()):
            future.cancel()

    curdoc().on_session_destroyed(session_closed)

    def _settle(doc, pid, name, result, on_done, exc):
        def apply():
            pending.pop(pid, None)
            pending_paths.pop(pid, None)
            pending_heavy.discard(pid)
            pending_futures.pop(pid, None)
            token = pending_tokens.pop(pid, None)
            if token is not None and token.cancelled:
                update_all(canvas)       # user cancelled: drop whatever the run returned
                return
            if exc is not None:
                logging.error("Failed to load %s: %s", name, exc)
            elif on_done is not None:
                try:
                    on_done(result)      # add_object -> on_change -> update_all
                except Exception:
                    logging.exception("Failed to add %s", name)
            update_all(canvas)           # clear the spinner (and show the new object)
        try:
            doc.add_next_tick_callback(apply)
        except Exception:
            pass  # session/document already gone

    def cancel_pending(pid):
        """Cancel a running background job from its row in the Processes panel. Drops the
        row now and cancels this workflow's dispatch subscriptions. Unshared subprocesses are killed;
        shared jobs continue for their other subscribers. Queued work is cancelled too.
        Either way the result (if any still arrives) is discarded in _settle."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        name = pending.pop(pid, None)
        pending_paths.pop(pid, None)
        pending_heavy.discard(pid)
        if name is None:
            return
        token = pending_tokens.get(pid)
        if token is not None:
            token.cancel()               # kill the dispatch process tree if it is running
        fut = pending_futures.get(pid)
        if fut is not None and fut.cancel():
            # still queued (never started): work()/_settle won't run, so clean up here
            pending_futures.pop(pid, None)
            pending_tokens.pop(pid, None)
        update_all(canvas)               # remove the spinner immediately
        notify(f"Cancelled {name}", level="info")

    def on_select(filepath):
        store.add_recent(filepath)
        name = Path(filepath).name
        load_async(name, lambda token=None: load_file(filepath),
                   lambda spec: canvas.add_object(spec, name=name,
                                                       source=("file", str(filepath)),
                                                       meta=getattr(spec, "meta", None)),
                   path=filepath)

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
        "select":       focus_spectrum,
        "color":        set_color,
        "profile":      ew.set_profile,
        "ew-snr":       ew.set_snr,
        "ew-cont":      ew.set_cont,
        "ew-toggle":    lambda v: tools.toggle("ew"),
        "ew-cancel":    lambda v: ew.deactivate(),
        "synth-toggle": lambda v: tools.toggle("synth"),
        "synth-cancel": lambda v: synth_panel.deactivate(),
        "synthesize":   lambda v: synth_panel.synthesize(),
        "annot-toggle": lambda v: ctx.annot.toggle_visible(),
        "linejump":     linenav.jump,
        "linesort":     linenav.set_sort,
        "lineclose":    lambda v: linenav.close(),
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
        sessions_bridge.value = sessions_html(store.list_sessions())

    def on_session_action(attr, old, new):
        if not new:
            return
        kind, _, rest = new.partition(" ")
        arg = rest.rsplit(" ", 1)[0]  # strip the trailing nonce
        if kind == "save":
            name = arg.strip() or time.strftime("Save %Y-%m-%d %H:%M")
            store.save_session(name, snapshot())
            render_sessions()
            notify(f"Saved “{name}”", level="success")
        elif kind == "load":
            snap = store.get_session(int(arg))
            if snap is None:
                notify("That save is no longer available.", level="danger")
                render_sessions()
            else:
                restore_workspace(snap)
        elif kind == "delete":
            store.delete_session(int(arg))
            render_sessions()
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
            synth_panel.delete_cache(arg)
        elif kind == "clear":
            synth_panel.clear_cache()
        elif kind == "refresh":
            synth_panel.refresh_cache()

    cache_actions.on_change("value", on_cache_action)

    # the "Load this run" button appears when the picker is inside a tumult run folder
    def on_picker_dir(folder):
        run_name_bridge.value = html.escape(Path(folder).name) if folder and is_run(folder) else ""
        toggle(run_ready_bridge, bool(folder and is_run(folder)))

    def load_current_run():
        folder = getattr(picker, "current_dir", None)
        if folder is None or not is_run(folder):
            return
        store.add_recent(folder)
        name = Path(folder).name

        def add(spec):
            color, dash = canvas.next_run_style()  # picked at add time -> no colour clashes
            canvas.add_object(spec, name=name, color=color, dash=dash,
                                   source=("run", str(folder)))

        load_async(name, lambda token=None: load_run(folder), add)

    # `busy_provider`: on_select hands the read to a worker thread and returns at once, so
    # the row's spinner is carried by the pending registry until the object lands
    picker = FilePicker(on_select, container_id="file-results", start_dir=start_dir,
                        on_dir=on_picker_dir, recent_provider=store.recent_files,
                        busy_provider=pending_paths.values)
    picker.layout.name = "file_picker"

    # a second FilePicker instance, embedded inline in the settings modal (no nested
    # modal) to change the navigator's master line list
    linelist_path_bridge = B.html("linelist-path")
    linelist_open_bridge = B.cls("#linelist-picker", "open")

    def on_linelist_select(filepath):
        strengths.set_master_path(filepath)
        if linenav.set_master_path(filepath):
            store.set_pref("master_linelist", str(filepath))
            ctx.annot.resync()
            linelist_path_bridge.value = html.escape(Path(filepath).name)
            pulse(linelist_open_bridge, "off")  # collapse the inline picker

    linelist_picker = FilePicker(on_linelist_select, container_id="linelist-results",
                                 start_dir=Path(strengths.master_path).parent)
    linelist_picker.layout.name = "linelist_picker"

    # every html control is a hidden model collected under one root: the bridges/actions
    # `B` gathered as they were created, plus each component's own form/bridge models
    bridges = column(*B.models,
                     *settings.models, *inspector.models, *linenav.models,
                     *linejump_form.models,
                     *synth_panel.models, *normtool.models, *filtertool.models,
                     *masktool.models,
                     *abundancetool.models, *abundcorrtool.models, *corr_atoms.models,
                     *ionbaltool.models)
    bridges.name = "bridges"

    doc = curdoc()
    doc.title = "TUMULT Studio"
    for root in (canvas.nav, canvas.main, filtertool.fig,
                 *abundancetool.figs, abundcorrtool.fig, *ionbaltool.figs, picker.layout,
                 *masktool.roots,
                 linelist_picker.layout, synth_panel.atoms.picker.layout,
                 synth_panel.grids.picker.layout, synth_panel.atmos_picker.layout,
                 abundancetool.setup.atmos_picker.layout,
                 abundcorrtool.setup.atmos_picker.layout,
                 corr_atoms.picker.layout, bridges):
        doc.add_root(root)

    ionbaltool.refresh()      # the balance modal's target selector

    inspector.bind(canvas.selected())
    sel = canvas.selected()
    active_color = sel.line.glyph.line_color if sel else None
    _auto = masktool.source("auto")     # None while masktool.AUTO_MASKS is off
    doc.template_variables.update(
        css=(HERE / "static" / "style.css").read_text(),
        main_empty=not canvas.objects,  # bake the empty-state class into the first paint
        objects_html=objects_panel_html(canvas),
        processes_html=processes_html(pending, pending_heavy),
        inspector_html=inspector_html(canvas),
        footer_html=footer_html(canvas),
        swatches_html=swatch_html(active_color),
        edit_html=inspector.edit_html(),
        style_html=inspector.style_form.html(),
        settings_html=settings.html(),
        synth_run_html=synth_panel.run_form.html(),
        synth_window_html=synth_panel.window_form.html(),
        synth_abund_html=synth_panel.abund_form.html(),
        synth_mode_html=synth_panel.mode_form.html(),
        synth_prefs_html=synth_panel.pref_forms["synthprefs"].html(),
        synth_engine_html=synth_panel.pref_forms["synthengine"].html(),
        grid_list_html=synth_panel.grids.list_html(),
        grid_results_name=synth_panel.grids.results_name,
        synth_grid_html=synth_panel.grid_ui_html(),
        synth_show_file=synth_panel.setup.show_file(),
        synth_file_results_name=synth_panel.atmos_picker.results_name,
        synth_cache_html=synth_panel.cache_html(),
        synth_warning=synth_panel.has_warning,
        synth_atoms_html=synth_panel.atoms.chips_html(),
        synth_is_nlte=synth_panel.is_nlte,
        atom_results_name=synth_panel.atoms.results_name,
        fit_fields_html=synth_panel.fit.html(),
        annot_label=ctx.annot.button_html(),
        linesettings_html=linejump_form.html(),
        linelist_path=html.escape(Path(strengths.master_path).name),
        linelist_results_name=linelist_picker.results_name,
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
        abund_atoms_ro_html=abund_setup.atoms_ro_html(),
        abund_file_results_name=abund_setup.atmos_picker.results_name,
        corr_results_html=abundcorrtool.results_html(),
        corr_table_html="",
        corr1_html=abundcorrtool.form1.html(),
        corr3_html=abundcorrtool.form3.html(),
        corr_params_html=corr_setup.run_form.html(),
        corr_grid_html=corr_setup.grid_ui_html(),
        corr_show_file=corr_setup.show_file(),
        corr_file_results_name=corr_setup.atmos_picker.results_name,
        corr_atoms_html=corr_atoms.chips_html(),
        corr_atom_results_name=corr_atoms.results_name,
        ionbal_results_html=ionbaltool.results_html(),
        ionbal_table_html="",
        ionbal_context_html="",
        ionbal1_html=ionbaltool.form1.html(),
        ionbal2_html=ionbaltool.form2.html(),
        ionbalsetup_html=ionbaltool.setup.form.html(),
        ionbaladv_html=ionbaltool.formadv.html(),
        form_script=FORM_SCRIPT + NUM_INPUT_SCRIPT,
        results_name=picker.results_name,
    )


