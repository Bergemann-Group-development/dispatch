"""TumultStudio2 — Bokeh server app entry point.

Run with:   bokeh serve --show TumultStudio2
"""
import html
import itertools
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

_nonce_counter = itertools.count()
_EXEC = ThreadPoolExecutor(max_workers=4)  # background loading (shared across sessions)
# synthesis runs one at a time: each build_stellar_grid spawns its own process pool
# sized to all cores, so concurrent jobs would oversubscribe the machine.
_SYNTH_EXEC = ThreadPoolExecutor(max_workers=1)


def _nonce():
    return next(_nonce_counter)

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))  # repo root → `import tumult` (during startup)

import numpy as np  # noqa: E402
from bokeh.io import curdoc  # noqa: E402
from bokeh.layouts import column  # noqa: E402

from spectrum import Spectrum, load_spectrum, load_file, is_run, load_run  # noqa: E402
from canvas import SpectrumCanvas, DASH_CSS, LINE  # noqa: E402
from context import AppContext  # noqa: E402
from controls import Field, Form, FORM_SCRIPT  # noqa: E402
from inspector import Inspector, swatch_html  # noqa: E402
from ew import EquivalentWidth  # noqa: E402
from filtertool import FilterTool  # noqa: E402
from normalisetool import NormaliseTool  # noqa: E402
from linenav import LineNavigator  # noqa: E402
from synthpanel import SynthPanel  # noqa: E402
from synthsetup import SynthSetup  # noqa: E402
from abundancetool import AbundanceTool  # noqa: E402
from abundcorrtool import AbundCorrTool  # noqa: E402
from sparamtool import SParamTool  # noqa: E402
import synth  # noqa: E402  (CancelToken/CancelledError for cancellable loads)
from filepicker import FilePicker  # noqa: E402
from toolmanager import ToolManager  # noqa: E402
from session import Workspace, ObjectSnap, FilterSnap  # noqa: E402
from htmlbridge import BridgeSet, pulse, toggle  # noqa: E402
from units import X_UNITS, Y_UNITS, to_range  # noqa: E402
from profiling import profiled  # noqa: E402
import store  # noqa: E402

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
    '<div class="object-empty">No spectra yet. Open one with the '
    '<i class="ti ti-folder-open"></i> button above, or click '
    '<b>Synth</b> to compute a new spectrum.</div>'
)


def objects_panel_html(canvas, pending, active_kind="spectrum"):
    """The Objects panel body: the loaded objects + any loading placeholders, or the
    empty-state hint when there is nothing to show."""
    rows = objects_html(canvas, active_kind) + pending_html(pending)
    return rows or OBJECTS_EMPTY_HINT


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


def pending_html(pending):
    """Spinner placeholders for objects still loading in the background, each with a
    cancel button (data-pid) that drops the run and discards its result."""
    return "".join(
        f'<div class="object-item pending" data-pid="{pid}">'
        '<span class="spinner"></span>'
        f'<span class="object-name" title="{html.escape(name)}">{html.escape(name)}</span>'
        f'<button class="object-btn danger" data-action="cancel-pending" data-uid="{pid}" '
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
    view_bridge = B.cls("#main", "lines-mode")  # swap nav overview / line table
    run_name_bridge = B.html("run-name")          # name shown on the "Load run" button
    run_ready_bridge = B.cls("#load-run", "ready")  # reveal it inside a run folder
    loading_bridge = B.cls("#footer", "loading")       # footer progress bar while loading
    empty_bridge = B.cls("#main", "empty")  # swap the empty canvas for the inline session loader
    synth_mode_bridge = B.cls(".sidebar-right", "synth-active")  # swap to the synth panel
    synth_btn_bridge = B.cls("#synth-btn", "active")             # toolbar button state
    synth_warn_bridge = B.cls("#settings-warning", "show")       # missing-files sign on the cog
    lines_panel_bridge = B.cls(".sidebar-right", "lines-active")  # line-nav inspector panel
    # (the abundance / correction wizards create and own their bridges via the context)

    # which left-sidebar list owns the single selection ("spectrum" | "result" | "corr");
    # the wizards drive it through ctx.focus_kind + the "focus-changed" event, and the
    # Objects list repaints to follow the highlight.
    ctx.focus_kind = "spectrum"
    ctx.subscribe("focus-changed",
                  lambda kind, **_: setattr(objects_bridge, "value",
                                            objects_panel_html(canvas, pending, kind)))

    def synth_on_mode(state):
        pulse(synth_mode_bridge, state)
        pulse(synth_btn_bridge, state)

    def synth_on_warn(has_warning):
        toggle(synth_warn_bridge, has_warning)

    inspector = linenav = abundancetool = abundcorrtool = sparamtool = None  # assigned after the canvas
    pending = ctx.pending  # id -> name of objects still loading in the background
    pending_futures = {}  # id -> Future, so a still-queued run can be cancelled before it starts
    pending_tokens = {}   # id -> synth.CancelToken, to hard-kill a running synthesis

    def run_busy(fn):
        """Run blocking work off the click callback so the footer progress bar paints
        first (the compute hogs the document thread, so a same-callback toggle never
        shows)."""
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
            objs.append(ObjectSnap(
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
            sparams=sparamtool.snapshot() if sparamtool is not None else {},
        ).to_dict()

    def update_all(canvas):
        objects_bridge.value = objects_panel_html(canvas, pending, ctx.focus_kind)
        inspector_bridge.value = inspector_html(canvas)
        footer_bridge.value = footer_html(canvas)
        synth_panel.refresh_targets()  # keep the Fit-tab target dropdown in sync
        toggle(loading_bridge, bool(pending))
        # show the inline session loader while the canvas holds nothing (and nothing
        # is on its way in)
        toggle(empty_bridge, not canvas.objects and not pending)
        if inspector is not None:
            inspector.bind(canvas.selected())
        ctx.emit("changed")  # EW / Filter / Normalise (and future tools) resync themselves
        if linenav is not None and linenav.active:
            linenav.render()  # keep the line list scoped to the active object's range
            # (only while the Lines view is open; it is rebuilt on re-entry otherwise)
        if abundancetool is not None:
            abundancetool.refresh()  # keep the step-1 reference selector in sync
        if sparamtool is not None:
            sparamtool.refresh()     # ditto for the Params wizard's target selector

    canvas = SpectrumCanvas(on_change=update_all)  # empty until a spectrum is opened
    ctx.canvas = canvas
    ctx.notify = notify
    ctx.run_busy = run_busy

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

        def run():
            try:
                tgt = target.spectrum.crop(xmin=lo, xmax=hi)
                best, _ = obj.spectrum.fit(tgt, fit=params, p0=dict(obj.conv))
            except Exception as exc:
                logging.exception("Edit-tab broadening fit failed")
                notify(f"Fit failed: {exc}")
                return
            obj.conv = {k: v for k, v in {**obj.conv, **best}.items() if v}
            canvas._refresh_object(obj)
            inspector.bind(obj)  # show the fitted values in the number fields
            update_all(canvas)
            notify("Broadening fit · "
                   + ", ".join(f"{k}={v:g}" for k, v in best.items()), level="success")

        run_busy(run)

    inspector = Inspector(canvas, on_fit=edit_fit, notify=notify)
    ew = EquivalentWidth(ctx)

    master = store.get_pref("master_linelist")
    if not master or not Path(master).is_file():
        master = MASTER_PATH
    linenav = LineNavigator(
        canvas, on_table=lambda h: setattr(linelist_bridge, "value", h), master_path=master)
    # the Lines view's own inspector panel: how far either side of a line to show on jump
    linejump_form = Form("linesettings", [
        Field("dist", "Distance from center [Å]", kind="number",
              value=linenav.jump_distance, min=0.01),
    ], on_change=lambda key, value: linenav.set_jump_distance(value))
    filtertool = FilterTool(ctx)
    normtool = NormaliseTool(ctx)
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
    abund_setup = SynthSetup(
        "abund", "abund_run", synth_panel.grids,
        atoms_provider=lambda: list(synth_panel.atoms.selected), notify=notify)
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

    # the stellar-parameter wizard: its own synthesis setup (whose run form is the solver's
    # initial guess), sharing the grid registry and the Synthesize tab's NLTE atoms.
    sparam_setup = SynthSetup(
        "sparam", "sparam_run", synth_panel.grids,
        atoms_provider=lambda: list(synth_panel.atoms.selected), notify=notify)
    sparamtool = SParamTool(
        ctx,
        setup=sparam_setup,
        list_observations=lambda: [(o.uid, o.name) for o in canvas.objects],
        get_observation=lambda uid: (canvas.effective(canvas._find(uid)) if canvas._find(uid)
                                     else None),
        master_path=master)

    # single-selection coordination between the Objects and the wizards' result lists
    def focus_spectrum(uid):
        if tools.active_name() in ("abundance", "abundcorr", "sparams"):
            tools.deactivate_active()  # selecting a spectrum leaves the wizard
        abundancetool.blur()          # leave the centre abund view, drop the result highlight
        abundcorrtool.blur()          # ditto for the correction view
        sparamtool.blur()             # ditto for the parameter view
        canvas.select_object(uid)     # -> on_change -> update_all repaints with spectrum focus

    # the abundance / parameter wizards' step-1 EWs track the reference's filter state
    ctx.subscribe("filter-state", lambda **_: abundancetool.on_filter_changed())
    ctx.subscribe("filter-state", lambda **_: sparamtool.on_filter_changed())

    # the mutually-exclusive sidebar tools: at most one owns the right panel at a time.
    # filter's enter is wrapped in run_busy so the progress bar paints during the FFT/fit.
    tools.register("ew", ew)
    tools.register("synth", synth_panel)
    tools.register("normalise", normtool)
    tools.register("filter", filtertool, runner=run_busy)
    tools.register("abundance", abundancetool)
    tools.register("abundcorr", abundcorrtool)
    tools.register("sparams", sparamtool)

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
                                       source=(obj.kind, obj.path or None), meta=meta)
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
        if sparamtool is not None:
            sparamtool.restore(ws.sparams)
        update_all(canvas)

    # --- background loading: placeholder + spinner now, real object when ready ----
    def load_async(name, loader, on_done, executor=_EXEC):
        pid = _nonce()
        pending[pid] = name
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

    # the wizards submit their synthesis work through the context (serialised on the
    # one-slot synth executor so concurrent grids don't oversubscribe the machine)
    ctx.submit = lambda name, loader, on_done: load_async(name, loader, on_done,
                                                          executor=_SYNTH_EXEC)

    def _settle(doc, pid, name, result, on_done, exc):
        def apply():
            pending.pop(pid, None)
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
        """Cancel a still-loading run from its sidebar spinner. Drops the spinner now and
        hard-kills the synthesis: the token tears down the dispatch process pool and its
        dispatch.x children; a run still queued behind another is cancelled before it starts.
        Either way the result (if any still arrives) is discarded in _settle."""
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            return
        name = pending.pop(pid, None)
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
                                                       meta=getattr(spec, "meta", None)))

    # client -> server channel for the sidebar hide/remove buttons ("action:uid:nonce")
    actions = B.action("object-actions")

    def set_color(value):
        obj = canvas.selected()
        if obj is not None:
            canvas.set_style(obj, color=value)

    def set_view(value):
        if value == "lines":
            # entering the Lines view: hand the right sidebar to the line-nav panel,
            # standing down whichever tool was occupying it
            tools.deactivate_active()
            pulse(view_bridge, "on")
            toggle(lines_panel_bridge, True)
            linenav.activate()
        else:
            linenav.deactivate()  # stop refreshing the hidden line table
            pulse(view_bridge, "off")
            toggle(lines_panel_bridge, False)
            canvas.reset_nav_view()  # the line nav zoomed it; restore the overview

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
        "view":         set_view,
        "linejump":     linenav.jump,
        "linesort":     linenav.set_sort,
        "species":      linenav.set_species,
        "loadrun":      lambda v: load_current_run(),
        "filter":       lambda v: tools.toggle("filter"),
        "filter-clear": lambda v: filtertool.clear(),
        "filtertoggle": lambda v: filtertool.set_visible(v == "on"),
        "normalise":    lambda v: tools.toggle("normalise"),
        "norm-init":    lambda v: run_busy(normtool.initialise),
        "norm-clear":   lambda v: normtool.clear(),
        "abundance":    lambda v: tools.toggle("abundance"),
        "abundcorr":    lambda v: tools.toggle("abundcorr"),
        "sparams":      lambda v: tools.toggle("sparams"),
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

    picker = FilePicker(on_select, container_id="file-results", start_dir=start_dir,
                        on_dir=on_picker_dir, recent_provider=store.recent_files)
    picker.layout.name = "file_picker"

    # a second FilePicker instance, embedded inline in the settings modal (no nested
    # modal) to change the navigator's master line list
    linelist_path_bridge = B.html("linelist-path")
    linelist_open_bridge = B.cls("#linelist-picker", "open")

    def on_linelist_select(filepath):
        if linenav.set_master_path(filepath):
            store.set_pref("master_linelist", str(filepath))
            linelist_path_bridge.value = html.escape(Path(filepath).name)
            pulse(linelist_open_bridge, "off")  # collapse the inline picker

    linelist_picker = FilePicker(on_linelist_select, container_id="linelist-results",
                                 start_dir=Path(linenav.master_path).parent)
    linelist_picker.layout.name = "linelist_picker"

    # every html control is a hidden model collected under one root: the bridges/actions
    # `B` gathered as they were created, plus each component's own form/bridge models
    bridges = column(*B.models,
                     *settings.models, *inspector.models, *linenav.models,
                     *linejump_form.models,
                     *synth_panel.models, *normtool.models, *filtertool.models,
                     *abundancetool.models, *abundcorrtool.models, *corr_atoms.models,
                     *sparamtool.models)
    bridges.name = "bridges"

    doc = curdoc()
    doc.title = "TUMULT Studio"
    for root in (canvas.nav, canvas.main, filtertool.fig, abundancetool.fig,
                 abundcorrtool.fig, sparamtool.fig, picker.layout,
                 linelist_picker.layout, synth_panel.atoms.picker.layout,
                 synth_panel.grids.picker.layout, synth_panel.atmos_picker.layout,
                 abundancetool.setup.atmos_picker.layout,
                 abundcorrtool.setup.atmos_picker.layout,
                 corr_atoms.picker.layout, sparamtool.setup.atmos_picker.layout, bridges):
        doc.add_root(root)

    inspector.bind(canvas.selected())
    sel = canvas.selected()
    active_color = sel.line.glyph.line_color if sel else None
    doc.template_variables.update(
        css=(HERE / "static" / "style.css").read_text(),
        main_empty=not canvas.objects,  # bake the empty-state class into the first paint
        objects_html=objects_panel_html(canvas, pending),
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
        linenav_html=linenav.species_html(),
        linesettings_html=linejump_form.html(),
        linelist_path=html.escape(Path(linenav.master_path).name),
        linelist_results_name=linelist_picker.results_name,
        sessions_html=sessions_html(store.list_sessions()),
        abund_results_html=abundancetool.results_html(),
        abund_table_html="",
        abund_ref_html=abundancetool.ref_html(),
        abund1_html=abundancetool.form1.html(),
        abund3_html=abundancetool.form3.html(),
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
        sparam_results_html=sparamtool.results_html(),
        sparam_table_html="",
        sparam_ref_html=sparamtool.ref_html(),
        sparam1_html=sparamtool.form1.html(),
        sparam2_html=sparamtool.form2.html(),
        sparam3_html=sparamtool.form3.html(),
        sparam_params_html=sparam_setup.run_form.html(),
        sparam_base_html=sparam_setup.base_form.html(),
        sparam_mode_html=sparam_setup.mode_form.html(),
        sparam_grid_html=sparam_setup.grid_ui_html(),
        sparam_show_file=sparam_setup.show_file(),
        sparam_atoms_ro_html=sparam_setup.atoms_ro_html(),
        sparam_file_results_name=sparam_setup.atmos_picker.results_name,
        form_script=FORM_SCRIPT,
        results_name=picker.results_name,
    )


build()
