"""TumultStudio document construction, callable by Bokeh and local harnesses.

Run with:   bokeh serve --show TumultStudio
"""
import html
import itertools
from copy import deepcopy
import logging
import os
import time
from pathlib import Path

import astropy.units as u
from bokeh.io import curdoc
from bokeh.layouts import column

from tumult.spectrum_io import Spectrum, load_file, is_run, load_run
from .spectrum.canvas import SpectrumCanvas, DASH_CSS, LINE
from .context import AppContext
from .ui.controls import Field, Form, FORM_SCRIPT, NUM_INPUT_SCRIPT
from .spectrum.inspector import Inspector, swatch_html
from .spectrum.filtertool import FilterTool
from .spectrum.normalisetool import NormaliseTool
from .masks.maskstore import MaskStore
from .masks.masktool import MaskTool
from .spectrum.linenav import LineNavigator
from .spectrum.lineannot import LineAnnotations
from tumult import refstars
from tumult import linestrength
from tumult.linestrength import LineStrengths
from .spectrum.rangetool import RangeTool
from .preferences.refstarspage import RefStarsPage
from .synthesis.synthmodal import SynthModal
from .synthesis.synthsetup import ATMOS_SCRIPT, SynthSetup
from .analysis.abundancetool import AbundanceTool
from .analysis.abundwizard import ANALYSIS_SCRIPT, PROFILES  # (shared EW-profile choices)
from .atoms.atomtool import AtomTool
from .analysis.abundcorrtool import AbundCorrTool
from .analysis.ionbalance import IonBalanceTool
from .workspace import naming
from .workspace.removal import Removal, register_workspace as register_removals
from .workspace.editmenu import EditMenu, Target as EditTarget
from tumult.synthesis import workflows as synth  # (execution_root of the configured prefs)
from .preferences.synthprefs import PREF_LABELS, load_prefs, _resolve_path  # (the configured workflow root)
from .ui.filepicker import FilePicker
from .ui.toolmanager import ToolManager
from .workspace.session import FOCUS_TOOL, RESULT_TOOLS
from .workspace.persistence import WorkspaceParts, WorkspacePersistence
from .ui.htmlbridge import BridgeSet, toggle
from .execution.delivery import DeliveryService
from .execution.processes import ProcessesPanel
from .execution.workflows import WorkflowService
from .ui.debounce import Trailing
from .ui.units import X_UNITS, Y_UNITS, to_range
from tumult import linestrength
from tumult.synthesis.runtime import (WORKFLOWS as _SYNTH_EXEC, HOP as _HOP_EXEC,
                                      get_runtime as _get_runtime)
from .workspace import store


HERE = Path(__file__).parent          # the TumultStudio package directory


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


def _load_atom(path):
    """Load a saved workspace's model atom (imported on demand: parsing an atom pulls in
    the whole atomic-data stack, which a document that opens none never needs)."""
    from tumult.atom import AtomClass
    return AtomClass(str(path))


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
    return (
        _field("Name", obj.name)   # renamed through the toolbar's Edit menu (editmenu.py)
        + _field("Range", f"{lo:.1f} – {hi:.1f} {canvas.xlabel}")
        + _field("Points", f"{len(obj.spectrum):,}")
        + _field("Native units", f"{obj.raw.xunit} / {obj.raw.yunit}")
        + _rerun_html(obj)
        + _meta_html(obj.meta)
    )


def _rerun_html(obj):
    """The Info tab's rerun action, for a spectrum that recorded its request."""
    if not (obj.meta or {}).get("rerun"):
        return ""
    return (f'<button class="text-btn object-rerun" data-uid="{html.escape(obj.uid)}" '
            'title="Open the synthesis modal with the request this spectrum was computed '
            'from">Rerun from this…</button>')


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
    # the job list, bottom left, and the daemon's live jobs it also shows
    panel = ProcessesPanel(ctx)
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

    def notify(message, level="warn", revert_ms=4000, undo=None):
        """``undo`` is a removal token (removal.py): the message then offers Undo."""
        lvl = level if level in _status_levels else "info"
        n = next(_status_token)
        _status_last["n"] = n
        button = ("" if undo is None else
                  f'<button class="status-action" data-undo="{html.escape(str(undo))}">'
                  'Undo</button>')
        # trailing comment keeps the pushed string unique, so an identical repeat
        # message still fires the bridge's change event and replays the flash
        status_bridge.value = (f'<span class="status-msg status-{lvl}">'
                               f'{html.escape(str(message))}{button}</span><!--{n}-->')

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
    edit_menu = None        # ditto: the toolbar's Edit menu, built once every list exists
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

    def focused():
        """``(kind, uid)`` of the object holding the single left-sidebar focus."""
        kind = ctx.focus_kind
        if kind == "spectrum":
            return kind, canvas.selected_uid
        if kind == "linemask":
            return kind, mask_store.current_uid
        if kind == "atom":
            return kind, atomtool.selected_uid
        return kind, wizards[FOCUS_TOOL[kind]].selected_uid

    def update_status():
        """Refresh job/loading surfaces without invalidating scientific tool state."""
        panel.repaint()
        footer_bridge.value = footer_html(canvas)
        active = bool(jobs.active)
        toggle(loading_bridge, active)
        toggle(empty_bridge, not canvas.objects
               and not any(job.delivers for job in jobs.active))
        if picker is not None:
            picker.refresh_busy()  # a load settled or was cancelled: stop its row spinning

    def update_all(canvas):
        objects_bridge.value = objects_panel_html(canvas, ctx.focus_kind)
        if edit_menu is not None:
            edit_menu.render()
        update_status()
        inspector_bridge.value = inspector_html(canvas)
        synth_modal.refresh_targets()  # keep the fit target dropdown in sync
        if inspector is not None:
            inspector.bind(canvas.selected())
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
        ctx.workspace_changed()      # the workspace changed: keep the autosave slot current

    canvas = SpectrumCanvas(on_change=update_all)  # empty until a spectrum is opened
    ctx.canvas = canvas
    ctx.notify = notify
    ctx.run_busy = run_busy
    ctx.busy = lambda on: toggle(loading_bridge, bool(on))
    # Background work: a Processes row now, the outcome delivered on the document. The
    # service installs ctx.submit / ctx.submit_io; `load_async` is its plain-pool entry
    # point, used by the loads and fits this module still composes.
    delivery = DeliveryService(ctx, panel, on_status=update_status)
    load_async = delivery.load_async
    ctx.execution = _get_runtime()

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
        # the target is read as measured (filtered when filtering is active); the fitted
        # object stays its own unconvolved spectrum, since the fit adjusts its convolution
        observed = canvas.target_spectrum(target.uid)
        lo, hi = to_range(r.start, r.end, canvas.xunit, observed.xunit)
        source, initial = obj.spectrum, dict(obj.conv)

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

    def duplicate_spectrum(obj):
        """Edit → Duplicate on a spectrum: ``obj`` again under a unique name, with its own
        copy of the data, provenance, continuum, convolution and style, selected so a
        broadening fit or a continuum edit changes the copy only."""
        g = obj.line.glyph
        copy = canvas.add_object(deepcopy(obj.raw), name=ctx.unique_name(obj.name, pending=False),
                                 color=g.line_color, dash=obj.dash, source=obj.source,
                                 meta=deepcopy(obj.meta), select=False)
        copy.norm_knots = deepcopy(obj.norm_knots)
        copy.already_normalised = obj.already_normalised
        copy.filtered = deepcopy(obj.filtered)
        copy.conv = dict(obj.conv)
        canvas.normalise_object(copy)
        canvas._refresh_object(copy)
        canvas.set_style(copy, width=g.line_width, alpha=g.line_alpha)
        canvas.select_object(copy.uid)
        notify(f"Duplicated as “{copy.name}”.", level="success")

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
    # the shared LineMask set (imported windows for fitting), owned centrally so
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
        get_observation=canvas.target_spectrum,
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
        o = canvas.add_object(spec, name=ctx.unique_name(name, pending=False),
                              color=color, dash=dash, source=("synth", None), meta=meta)
        if conv:
            o.conv = dict(conv)
            canvas._refresh_object(o)  # update_all (in _settle) then rebinds the Edit tab
        synth_modal.refresh_cache()  # the new run is now in the cache listing

    # This document's attachment to the durable workflows every dispatch calculation is:
    # the service installs ctx.submit_operation and the recovery actions, so every
    # synthesis and every fit built below submits through the context.
    workflows = WorkflowService(ctx, delivery, panel, add_spectrum=add_synth,
                                workflow_roots=workflow_roots, on_status=update_status)
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
    from .synthesis.atompicker import AtomPicker
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
        get_observation=canvas.target_spectrum,
        master_path=master)

    # the abundance-correction tool: its own synthesis setup (target/NLTE atmosphere) in its
    # modal, sharing only the app-global grid registry (one grid flagged as the 1D LTE
    # reference). The model atom is the atom object the modal was opened from.
    corr_setup = SynthSetup(
        "corr", "abundcorr_run", synth_modal.grids,
        atoms_provider=lambda: ([abundcorrtool._pending_atom.path]
                                if abundcorrtool._pending_atom is not None else []),
        notify=notify)

    abundcorrtool = AbundCorrTool(
        ctx,
        setup=corr_setup,
        grids=synth_modal.grids)

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
        get_observation=canvas.target_spectrum,
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

    # --- naming, removal and the Edit menu, across every object list ------
    # the three result lists, by the name of their tool (the session's result sections)
    wizards = {"abundance": abundancetool, "abundcorr": abundcorrtool, "ionbal": ionbaltool}
    assert tuple(wizards) == RESULT_TOOLS

    def unique_name(name, pending=True, own=None):
        """``name`` made unique among the spectrum and result objects (:func:`naming.unique`);
        ``pending`` also avoids the names of this document's unfinished jobs, so two
        submissions of the same calculation deliver distinguishable objects. ``own`` is the
        current name of an object being renamed, which does not count as taken."""
        taken = {o.name for o in canvas.objects}
        taken.update(r.name for w in wizards.values() for r in w.results)
        if pending:
            taken.update(job.name for job in jobs.active)
        taken.discard(own)
        return naming.unique(name, taken)

    ctx.unique_name = unique_name

    # every sidebar remove button goes through this: it asks when the removal would cost
    # something and offers Undo in the status message (removal.py)
    ctx.removal = Removal(ctx, notify)
    register_removals(ctx.removal, canvas, masktool, wizards, atomtool)

    def rename_spectrum(obj, name):
        canvas.rename_object(obj.uid, unique_name(name, own=obj.name))

    edit_menu = ctx.edit_menu = EditMenu(ctx, focused)
    edit_menu.register("spectrum", EditTarget(
        find=canvas._find, removal="spectrum",
        rename=rename_spectrum, duplicate=duplicate_spectrum))

    def duplicate_set(mset):
        copy = mask_store.duplicate(mset.uid)
        if copy is not None:
            masktool.select_set(copy.uid)
            notify(f"Duplicated as “{copy.name}”.", level="success")

    edit_menu.register("linemask", EditTarget(
        find=ctx.removal.kinds["linemask"].find, removal="linemask",
        rename=lambda mset, name: mask_store.rename(mset.uid, name),
        duplicate=duplicate_set))
    for kind, tool_name in FOCUS_TOOL.items():
        wizard = wizards[tool_name]
        edit_menu.register(kind, EditTarget(
            find=wizard._find, removal=tool_name,
            rename=lambda res, name, w=wizard: w.rename_result(res.uid, name),
            duplicate=lambda res, w=wizard: w.duplicate_result(res.uid)))
    edit_menu.register("atom", EditTarget(
        find=atomtool.find, removal="atom",
        fixed_note="A model atom is identified by its file"))
    edit_menu.render()

    # --- session save / restore -------------------------------------------
    # The workspace owner: the snapshot, the validate -> prepare -> commit flow with its
    # rollback, the skipped-file report and the autosave slot. It installs
    # ctx.restore_workspace (local harnesses load slots through it) and
    # ctx.workspace_changed, and owns the session-slot action channel and its dialogs.
    persist = WorkspacePersistence(
        ctx,
        WorkspaceParts(atoms=atomtool, mask_view=masktool, results=wizards,
                       filter=filtertool, settings=settings, focused=focused,
                       focus_spectrum=focus_spectrum,
                       refresh=lambda: update_all(canvas)),
        delivery, workflows, load_spectrum=_load_obj, load_atom=_load_atom)

    def use_run_departures(run_id):
        """Hand a departure solve to the synthesis modal's Saved departures field."""
        entry = workflows.external_run(run_id)
        if entry is None:
            return
        synth_modal.setup.mode_form.set_values(mode="NLTE", departures=entry.folder)
        synth_modal.setup._save_mode("mode", "NLTE")
        synth_modal.setup._save_mode("departures", entry.folder)
        notify(f"Synthesis modal set to reuse departures from {entry.name}", level="info")

    def session_closed(session_context):
        delivery.cancel_document_work()   # loads, fits, analysis jobs; not the detached
        # Cancel a queued older autosave before the synchronous durability boundary. If it
        # is already writing, the store lock orders this newer snapshot after it.
        persist.autosave(now=True)   # the document's last state, for Resume last session
        ctx.closed = True
        panel.live.close()

    ctx.doc.on_session_destroyed(session_closed)

    def on_select(filepath):
        remember_recent(filepath)
        name = Path(filepath).name
        def focus_state():
            return (ctx.focus_kind, canvas.selected_uid, tools.active_name(),
                    canvas.main.x_range.start, canvas.main.x_range.end)
        initial = focus_state()
        def add_loaded(spec):
            select = focus_state() == initial
            canvas.add_object(spec, name=unique_name(name, pending=False),
                              source=("file", str(filepath)),
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
        "remove":       lambda uid: ctx.removal.request("spectrum", uid),
        "cancel-pending": delivery.cancel_pending,
        "dismiss-job":  delivery.dismiss_job,
        "workflow-discover": lambda v: workflows.discover(),
        "workflow-attach": workflows.attach,
        "workflow-import": workflows.import_result,
        "workflow-forget": workflows.forget,
        "workflow-cancel": workflows.cancel_workflow,
        "job-details":  workflows.show_details,
        "folder-open":  workflows.open_folder,
        "run-cancel":   workflows.cancel_run,
        "run-departures": use_run_departures,
        "select":       focus_spectrum,
        "color":        set_color,
        "profile":      lambda v: range_tool.ew_setting(range_tool.ew.set_profile, v),
        "ew-snr":       lambda v: range_tool.ew_setting(range_tool.ew.set_snr, v),
        "ew-cont":      lambda v: range_tool.ew_setting(range_tool.ew.set_cont, v),
        "range-ew":     lambda v: range_tool.toggle_ew(),
        "mask-toggle":  lambda v: tools.toggle("range"),
        "range-cancel": lambda v: range_tool.cancel(),
        "escape":       lambda v: tools.cancel_active(),
        "range-new-mask": lambda v: range_tool.new_mask(),
        "range-synth":  lambda v: range_tool.synthesise(),
        "range-fit":    lambda v: range_tool.fit(),
        "range-add":    lambda v: range_tool.add_window(),
        "synth-new":    lambda v: synth_modal.open_range(),
        "rerun":        lambda v: synth_modal.open_rerun(canvas._find(v)),
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
        "norm-cancel":  lambda v: normtool.cancel(),
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
            canvas.add_object(spec, name=unique_name(name, pending=False),
                              color=color, dash=dash, source=("run", str(folder)))

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

    from .synthesis.draft import CalculationDraft
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
    # each add_root recomputes the whole reachable-model set from scratch (bokeh
    # Document.models.recompute); with ~25 roots that is quadratic in the final model
    # count. `models.freeze()` defers that recomputation to a single pass after every
    # root is in place, which is the difference between milliseconds and seconds here.
    with doc.models.freeze():
        for root in (canvas.nav, canvas.main, filtertool.fig,
                     *abundancetool.figs, abundancetool.cog.fig, abundancetool.chi2.fig,
                     *abundcorrtool.figs, abundcorrtool.cog.fig,
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

    panel.repaint()           # the daemon's jobs are rows from the first paint on

    inspector.bind(canvas.selected())
    sel = canvas.selected()
    active_color = sel.line.glyph.line_color if sel else None
    _auto = masktool.source("auto")     # None while masktool.AUTO_MASKS is off
    doc.template_variables.update(
        css=(HERE / "static" / "style.css").read_text(),
        main_empty=not canvas.objects,  # bake the empty-state class into the first paint
        objects_html=objects_panel_html(canvas),
        processes_html=panel.body_html,
        processes_usage_html=panel.usage_html,
        inspector_html=inspector_html(canvas),
        footer_html=footer_html(canvas),
        swatches_html=swatch_html(active_color),
        edit_html=inspector.edit_html(),
        edit_menu_html=edit_menu.html(),
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
        sessions_html=persist.slots_html(),
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
