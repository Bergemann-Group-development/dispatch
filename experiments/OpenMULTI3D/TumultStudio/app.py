"""TumultStudio document construction, callable by Bokeh and local harnesses.

Run with:   python -m TumultStudio   (or bokeh serve --show TumultStudio, without the landing page)

:func:`build` composes one document: the :class:`~.context.AppContext` with its services
(delivery, workflows, persistence), the workspace and the canvas over it, the tools and
views, the action channels that route the page's gestures to their owners, and the Bokeh
roots. Behaviour lives in those owners; what is here is construction, the connections
between them and the template variables of the first paint.
"""
import html
import logging
from pathlib import Path

from bokeh.io import curdoc
from bokeh.layouts import column

from tumult.synthesis import workflows as synth  # (execution_root of the configured prefs)
from tumult.synthesis.runtime import get_runtime

from .analysis import results as analysis_results
from .analysis.abundancetool import AbundanceTool
from .analysis.abundcorrtool import AbundCorrTool
from .analysis.ionbalance import IonBalanceTool
from .analysis.resultview import ANALYSIS_SCRIPT
from .atoms.atomtool import AtomTool
from .atoms.departuretool import DepartureTool
from .context import AppContext
from .execution.delivery import DeliveryService
from .execution.processes import ProcessesPanel
from .execution.workflows import WorkflowService
from .masks.maskstore import MaskStore
from .masks.masktool import MaskTool
from .preferences.page import PreferencesPage
from .preferences.refstarspage import RefStarsPage
from .preferences.synthprefs import PREF_LABELS, REPO_DIR, load_prefs, _resolve_path
from .spectrum.canvas import SpectrumCanvas, settings_form
from .spectrum.ew import FIT_PROFILES  # (the EW profiles an abundance may be measured with)
from .spectrum.filtertool import FilterTool
from .spectrum.inspector import Inspector, swatch_html
from .spectrum.normalisetool import NormaliseTool
from .spectrum.rangetool import RangeTool
from .spectrum.strengths import LineStrengthService
from .synthesis.atompicker import AtomPicker
from .synthesis.draft import CalculationDraft
from .synthesis.synthmodal import SynthModal
from .synthesis.synthsetup import ATMOS_SCRIPT, SynthSetup
from .ui.controls import FORM_SCRIPT, NUM_INPUT_SCRIPT
from .ui.htmlbridge import BridgeSet, StatusLine, toggle
from .ui.toolmanager import ToolManager
from .workspace import store
from .workspace.editmenu import EditMenu, Target as EditTarget
from .workspace.objects import WorkspaceObjects
from .workspace.opening import SpectrumOpener
from .workspace.persistence import WorkspacePersistence
from .workspace.removal import Removal, register_workspace as register_removals
from .workspace.session import FOCUS_TOOL


HERE = Path(__file__).parent          # the TumultStudio package directory

# the page's style sheets, one per owner, inlined into the template in this order: base
# first (tokens, reset and the shared controls), then the shell, then the dialogs. The
# order is the cascade, so a feature sheet may refine a base rule of equal specificity.
CSS_SHEETS = ("base", "workspace", "preferences", "synthesis", "analysis")


def page_css() -> str:
    """The inlined style text: the sheets in cascade order."""
    return "\n".join((HERE / "static" / "css" / f"{name}.css").read_text()
                     for name in CSS_SHEETS)


# the page's scripts, loaded from the app's static route in this order (see index.html)
JS_SCRIPTS = ("bridge", "workspace", "preferences", "synthesis", "analysis")


def js_version() -> str:
    """A query string that changes whenever a script does, so a browser refetches it.

    The static route sends no Cache-Control, so a browser may reuse a cached script on a
    plain navigation (heuristic freshness) while the inlined CSS and the templates come
    fresh from a restarted server; the two then disagree until a hard reload.
    """
    return str(max(int((HERE / "static" / "js" / f"{name}.js").stat().st_mtime)
                   for name in JS_SCRIPTS))


logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(name)s: %(message)s")
# Job failures are logged from the worker thread, traceback included, whether or not the
# document that admitted them still exists; the user-facing report goes through the
# owning document's Processes row and footer alone. No root-logger handler is installed
# per document, so one session's errors never reach another session's footer.


def workflow_roots():
    """Where durable workflow records are looked for: the synthesis output root the
    preferences configure (``<dispatch_cwd>/<synth_output_dir>``), which is where
    ``tumult.synthesis.submit`` writes ``_workflows/``. Discovery never leaves it."""
    prefs = load_prefs()
    return [synth.execution_root(prefs) / (prefs.get("synth_output_dir") or "synth")]


def build():
    # every hidden bridge / action channel made below is collected here, so the document
    # root column at the end is just `*B.models` instead of a hand-listed enumeration.
    # The context carries it (plus the ToolManager and, later, canvas/notify/run_busy)
    # into the tools, which create and own their bridges themselves.
    B = BridgeSet()
    tools = ToolManager()
    ctx = AppContext(curdoc(), B, tools)
    workspace = ctx.workspace = WorkspaceObjects(ctx)   # the live objects and the focus

    objects_bridge = B.html("objects")
    # the job list, bottom left, and the daemon's live jobs it also shows
    panel = ProcessesPanel(ctx)
    inspector_bridge = B.html("inspector")
    footer_bridge = B.html("footer-info")
    # transient info messages (left of footer); named for the same reason
    status = StatusLine(B.html("footer-status", name="footer-status"), ctx.doc)
    ctx.notify = notify = status.notify
    loading_bridge = B.cls("#footer", "loading")       # footer progress bar while loading
    empty_bridge = B.cls("#main", "empty")  # swap the empty canvas for the inline session loader
    synth_warn_bridge = B.cls("#settings-warning", "show")       # missing-files sign on the cog

    # the Objects list repaints to follow the single left-sidebar highlight
    ctx.subscribe("focus-changed",
                  lambda kind, **_: setattr(objects_bridge, "value",
                                            workspace.objects_panel_html(kind)))

    inspector = lines = abundancetool = ionbaltool = None  # assigned after the canvas
    edit_menu = None        # ditto: the toolbar's Edit menu, built once every list exists
    opener = None           # ditto: the spectrum file picker, whose rows spin while loading
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

    # --- the document's repaints: who follows a change of the objects or of the jobs ---
    def update_status():
        """Refresh job/loading surfaces without invalidating scientific tool state."""
        panel.repaint()
        footer_bridge.value = workspace.footer_html()
        active = bool(jobs.active)
        toggle(loading_bridge, active)
        toggle(empty_bridge, not workspace.spectra.objects
               and not any(job.delivers for job in jobs.active))
        if opener is not None:
            opener.picker.refresh_busy()  # a load settled or was cancelled: stop its row spinning

    def update_all(canvas=None):
        objects_bridge.value = workspace.objects_panel_html()
        if edit_menu is not None:
            edit_menu.render()
        update_status()
        inspector_bridge.value = inspector.info_html()
        synth_modal.refresh_targets()  # keep the fit target dropdown in sync
        inspector.bind(workspace.spectra.selected())
        ctx.emit("changed")  # EW / Filter / Normalise (and future tools) resync themselves
        if lines is not None:
            lines.resync()       # a different object may be a different star (linestrength)
        if abundancetool is not None:
            abundancetool.refresh()  # keep the abundance modal's reference selector in sync
        if ionbaltool is not None:
            ionbaltool.refresh()     # ditto for the balance modal's target selector
        ctx.workspace_changed()      # the workspace changed: keep the autosave slot current

    canvas = SpectrumCanvas(workspace.spectra, on_change=update_all)  # empty until a spectrum is opened
    ctx.canvas = canvas
    workspace.on_change = update_all
    ctx.run_busy = run_busy
    ctx.busy = lambda on: toggle(loading_bridge, bool(on))
    # Background work: a Processes row now, the outcome delivered on the document. The
    # service installs ctx.submit / ctx.submit_io; `load_async` is its plain-pool entry
    # point, used by the loads composed below.
    delivery = DeliveryService(ctx, panel, on_status=update_status)
    load_async = delivery.load_async
    ctx.execution = get_runtime()

    # the Select range tool: the one range gesture, whose actions reach the line-mask view and the
    # synthesis modal built further down (the lambdas resolve them at click time)
    inspector = Inspector(ctx, canvas,
                          on_select_range=lambda lo, hi: range_tool.select_range(lo, hi))
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
    # Preferences (sheets, pickers, the app-global grid registry, the run cache) and the
    # shared NLTE model-atom selection are built here, before the calculation dialogs that
    # consume them; on_inputs_change / on_grids_change bind once their listeners exist.
    prefs_page = PreferencesPage(ctx, on_warn=lambda on: toggle(synth_warn_bridge, on))
    _cwd = Path(str(prefs.get("dispatch_cwd") or REPO_DIR)).expanduser()
    _synth_atoms_dir = _cwd / "input" / "atoms"     # the synthesis modal's atom selector
    synth_atoms = AtomPicker(
        start_dir=_synth_atoms_dir if _synth_atoms_dir.is_dir() else _cwd)
    # the line strengths: the ranking the annotations and the line-list panel read
    lines = LineStrengthService(ctx, master, load_async)
    ctx.strengths = lines.strengths
    ctx.annot = lines.annot

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
        list_observations=workspace.spectra.listed,
        get_observation=canvas.target_spectrum,
        master_path=master, profiles=FIT_PROFILES,
        # "+ New mask" draws its window with the Select range tool's gesture
        place_window=range_tool.begin_add,
        # the mask view's analysis actions: hand the open set to the abundance / balance
        # tools and the synthesis modal, which are built further down (the lambdas resolve
        # them at click time)
        on_fit_balance=lambda mset: ionbaltool.open_for(mset),
        on_measure_abund=lambda mset: abundancetool.open_for(mset),
        on_synth=lambda mset: synth_modal.open_masks(mset, "synth"),
        on_fit=lambda mset: synth_modal.open_masks(mset, "fit"))
    settings = settings_form(canvas)

    def add_synth(spec, name, meta=None, conv=None):
        workspace.add_synthetic(spec, name, meta=meta, conv=conv)
        prefs_page.refresh_cache()  # the new run is now in the cache listing

    # This document's attachment to the durable workflows every dispatch calculation is:
    # the service installs ctx.submit_operation and the recovery actions, so every
    # synthesis and every fit built below submits through the context.
    workflows = WorkflowService(ctx, delivery, panel, add_spectrum=add_synth,
                                workflow_roots=workflow_roots, on_status=update_status)
    ctx.workflows = workflows          # local harnesses follow jobs through it
    synth_modal = SynthModal(ctx, add=add_synth, grids=prefs_page.grids, atoms=synth_atoms)
    ctx.synth_modal = synth_modal      # local harnesses drive the modal through it
    lines.use_marcs_grid(prefs_page.grids.grids)

    # the Abundance tool owns its synthesis setup (its own Setup sub-tab), sharing only the
    # app-global grid registry and the synthesis modal's NLTE model-atom selection; it inherits
    # the modal's configuration on first open and is independent thereafter.
    # its fields live in the abundance modal, so every edit restates the setup in the object
    # panel, where the Run synthesis button stands
    _abund_note = lambda *_: (abundancetool is not None and abundancetool.push_setup_note())
    abund_atoms = AtomPicker(Path(load_prefs()["dispatch_cwd"]) / "input" / "atoms", ns="abund-atoms", section="abund_run",
        picker_sel="#abund-atom-picker", section_sel="#abund-atoms",
        actions_name="abund-atom-actions", results_id="abund-atom-results",
        layout_name="abund_atom_picker", del_class="abund-atom-del")
    abund_setup = SynthSetup(
        "abund", "abund_run", prefs_page.grids,
        atoms_provider=lambda: list(abund_atoms.selected), notify=notify,
        on_run_change=_abund_note, on_base_change=_abund_note,
        on_mode_change=lambda on: (abund_atoms.show(on), _abund_note()), on_atmos_change=_abund_note)
    abundancetool = AbundanceTool(
        ctx,
        setup=abund_setup,
        list_observations=workspace.spectra.listed,
        get_observation=canvas.target_spectrum,
        master_path=master)

    # the abundance-correction tool: its own synthesis setup (target/NLTE atmosphere) in its
    # modal, sharing only the app-global grid registry (one grid flagged as the 1D LTE
    # reference). The model atom is the atom object the modal was opened from.
    corr_setup = SynthSetup(
        "corr", "abundcorr_run", prefs_page.grids,
        atoms_provider=lambda: ([abundcorrtool._pending_atom.path]
                                if abundcorrtool._pending_atom is not None else []),
        notify=notify)

    abundcorrtool = AbundCorrTool(
        ctx,
        setup=corr_setup,
        grids=prefs_page.grids)

    # the departure-coefficient tool: its own atmosphere setup in its modal, the model atom
    # fixed by the object the modal was opened from. A solved run is handed to the synthesis
    # modal's Saved departures field, which is what makes the populations reusable.
    def use_saved_departures(folder):
        synth_modal.setup.use_departures(folder)
        notify(f"Synthesis modal set to reuse the departures in {Path(folder).name}")

    dep_setup = SynthSetup(
        "dep", "dep_run", prefs_page.grids,
        atoms_provider=lambda: ([departuretool._pending_atom.path]
                                if departuretool._pending_atom is not None else []),
        notify=notify, allow_interpolation=False)
    departuretool = DepartureTool(ctx, setup=dep_setup, load_async=load_async,
                                  use_departures=use_saved_departures)

    # the model-atom objects: Open → Model atom loads one; its view is the two
    # tables, its panel the entry point of the departure and correction modals, and both
    # kinds of result are listed under it
    atomtool = AtomTool(ctx, corr=abundcorrtool, dep=departuretool)
    abundcorrtool.on_results = atomtool.render_objects
    abundcorrtool._known_atom_paths = lambda: [a.path for a in atomtool.atoms]
    departuretool.on_results = atomtool.render_objects
    departuretool._known_atom_paths = lambda: [a.path for a in atomtool.atoms]
    departuretool.atom_object = atomtool.by_path

    # The ionisation-balance modal owns its grid selection and always interpolates it.  Only
    # the app-global registered-grid list and synthesis preferences are shared.
    ionbaltool = IonBalanceTool(
        ctx,
        grids=prefs_page.grids,
        list_observations=workspace.spectra.listed,
        get_observation=canvas.target_spectrum,
        master_path=master)

    def on_grids_change():
        """The grid registry changed (Preferences → Atmosphere grids): every setup re-lists
        its grid, and the strength ranking follows the MARCS grid. The page has already
        refreshed its own file warnings."""
        synth_modal.setup.refresh_grids()
        abund_setup.refresh_grids()
        corr_setup.refresh_grids()
        dep_setup.refresh_grids()
        lines.use_marcs_grid(prefs_page.grids.grids)
        lines.warm()

    prefs_page.on_grids_change = on_grids_change

    def on_refstars_change():
        """The reference-star catalogue was edited: repaint every dropdown listing it."""
        for setup in (synth_modal.setup, abund_setup, corr_setup, dep_setup):
            setup.refresh_refstars()
        ionbaltool.refresh_refstars()
        lines.refresh_refstars()

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
    tools.register("departures", departuretool)

    # --- naming, removal and the Edit menu, across every object list ------
    # the three result lists, by the name of their tool (the session's result sections)
    wizards = {"abundance": abundancetool, "abundcorr": abundcorrtool, "ionbal": ionbaltool}
    workspace.bind(canvas, masks=mask_store, mask_view=masktool, atoms=atomtool,
                   results=wizards, departures=departuretool)

    # every sidebar remove button goes through this: it asks when the removal would cost
    # something and offers Undo in the status message (removal.py)
    ctx.removal = Removal(ctx, notify)
    register_removals(ctx.removal, workspace)

    edit_menu = ctx.edit_menu = EditMenu(ctx, workspace.focused)
    edit_menu.register("spectrum", EditTarget(
        find=workspace.spectra.find, removal="spectrum",
        rename=workspace.rename_spectrum, duplicate=workspace.duplicate_spectrum))
    edit_menu.register("linemask", EditTarget(
        find=ctx.removal.kinds["linemask"].find, removal="linemask",
        rename=workspace.rename_set, duplicate=workspace.duplicate_set))
    for kind, tool_name in FOCUS_TOOL.items():
        wizard = wizards[tool_name]
        edit_menu.register(kind, EditTarget(
            find=wizard.find, removal=tool_name,
            rename=lambda res, name, w=wizard: w.rename_result(res.uid, name),
            duplicate=lambda res, w=wizard: w.duplicate_result(res.uid)))
    edit_menu.register("atom", EditTarget(
        find=atomtool.find, removal="atom",
        fixed_note="A model atom is identified by its file"))
    edit_menu.register("departures", EditTarget(
        find=departuretool.find, removal="departures",
        rename=lambda obj, name: departuretool.rename(obj.uid, name),
        fixed_note="A departure solve is its run on disk"))
    edit_menu.render()

    # --- session save / restore -------------------------------------------
    # The workspace owner: the snapshot, the validate -> prepare -> commit flow with its
    # rollback, the skipped-file report and the autosave slot. It installs
    # ctx.restore_workspace (local harnesses load slots through it) and
    # ctx.workspace_changed, and owns the session-slot action channel and its dialogs.
    persist = WorkspacePersistence(
        ctx, workspace, delivery, workflows,
        filter=filtertool, settings=settings, refresh=update_all,
        deserialise_results=analysis_results.deserialise,
        relink_results=analysis_results.relink)

    def use_run_departures(run_id):
        """Hand a departure solve to the synthesis modal's Saved departures field."""
        entry = workflows.external_run(run_id)
        if entry is None:
            return
        synth_modal.setup.use_departures(entry.folder)
        notify(f"Synthesis modal set to reuse departures from {entry.name}", level="info")

    def session_closed(session_context):
        delivery.cancel_document_work()   # loads, fits, analysis jobs; not the detached
        # Cancel a queued older autosave before the synchronous durability boundary. If it
        # is already writing, the store lock orders this newer snapshot after it.
        persist.autosave(now=True)   # the document's last state, for Resume last session
        ctx.closed = True
        panel.live.close()

    ctx.doc.on_session_destroyed(session_closed)

    # client -> server channel for the sidebar hide/remove buttons ("action:uid:nonce")
    actions = B.action("object-actions")

    def set_color(value):
        obj = workspace.spectra.selected()
        if obj is not None:
            canvas.set_style(obj, color=value)

    # clicking a line-position indicator opens the line list in the bottom half of the
    # workspace, scrolled to that transition's row (lineannot -> "line-picked" -> linenav)
    linenav = lines.linenav
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
        "select":       workspace.focus_spectrum,
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
        "rerun":        lambda v: synth_modal.open_rerun(workspace.spectra.find(v)),
        "synth-submit": lambda v: synth_modal.submit(),
        "annot-toggle": lambda v: ctx.annot.toggle_visible(),
        "linejump":     linenav.jump,
        "linesort":     linenav.set_sort,
        "lineclose":    lambda v: linenav.close(),
        "lineopen":     lambda v: linenav.close() if linenav.active else linenav.open_at(),
        "loadrun":      lambda v: opener.open_run(),
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
        "departures":   lambda v: tools.toggle("departures"),
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

    # Open → Spectrum (the file picker and the "Load this run" button) and Open → Model atom
    opener = SpectrumOpener(ctx, workspace, load_async)
    atom_picker = atomtool.make_picker(load_async)

    def sync_linelist():
        prefs = load_prefs()
        path = str(_resolve_path(prefs, prefs["linelist"]))
        abundancetool.master_path = path
        ionbaltool.master_path = path
        auto = masktool.source("auto")
        if auto is not None:
            auto.linelist_path = path
        lines.set_master_path(path)

    prefs_page.on_inputs_change = sync_linelist
    if store.get_pref("master_linelist") is not None:
        with store.transaction() as data:
            data["preferences"].pop("master_linelist", None)
        notify("Navigation and analysis now use the synthesis Line list in Preferences.")

    # each calculation dialog's edits are a draft until submitted; closing the dialog
    # discards it (the "dialog-actions" channel names the dialog)
    synth_modal.draft = CalculationDraft(setup=synth_modal.setup, atoms=synth_atoms,
        owner=synth_modal.fit, attributes=("vals", "flags", "abund_rows", "target"))
    abundancetool.draft = CalculationDraft(forms=[("abundance", abundancetool.form1)],
        setup=abundancetool.setup, atoms=abund_atoms, owner=abundancetool, attributes=("ref_uid",))
    abundcorrtool.draft = CalculationDraft(forms=[("abundcorr", abundcorrtool.form1)],
        setup=abundcorrtool.setup)
    departuretool.draft = CalculationDraft(forms=[("departures", departuretool.form)],
        setup=dep_setup)
    ionbaltool.draft = CalculationDraft(forms=[("ionbal", form) for form in
        (ionbaltool.form1, ionbaltool.form2, ionbaltool.formadv, ionbaltool.setup.form)])
    drafts = {"synth-modal": synth_modal.draft, "abund-modal": abundancetool.draft,
              "corr-modal": abundcorrtool.draft, "dep-modal": departuretool.draft,
              "ionbal-modal": ionbaltool.draft}
    dialog_actions = B.action("dialog-actions")
    def close_draft(attr, old, new):
        draft = drafts.get(new.split(" ")[0])
        if draft is not None:
            draft.discard()
    dialog_actions.on_change("value", close_draft)

    # every html control is a hidden model collected under one root: the bridges/actions
    # `B` gathered as they were created, plus each component's own form/bridge models
    bridges = column(*B.models,
                     *settings.models, *inspector.models, *lines.models,
                     *synth_atoms.models, *abund_atoms.models,
                     *refstars_page.models, *prefs_page.models,
                     *synth_modal.models, *range_tool.models, *normtool.models,
                     *filtertool.models, *masktool.models,
                     *abundancetool.models, *abundcorrtool.models,
                     *departuretool.models, *ionbaltool.models)
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
                     *abundcorrtool.figs, abundcorrtool.cog.fig, *departuretool.figs,
                     *ionbaltool.figs, opener.picker.layout,
                     *masktool.roots,
                     synth_atoms.picker.layout, abund_atoms.picker.layout,
                     *(p.layout for p in prefs_page.pref_pickers.values()),
                     prefs_page.grids.picker.layout, synth_modal.atmos_picker.layout,
                     abundancetool.setup.atmos_picker.layout,
                     abundcorrtool.setup.atmos_picker.layout,
                     dep_setup.atmos_picker.layout,
                     atom_picker.layout, bridges):
            doc.add_root(root)

    ionbaltool.refresh()      # the balance modal's target selector

    panel.repaint()           # the daemon's jobs are rows from the first paint on

    inspector.bind(workspace.spectra.selected())
    sel = workspace.spectra.selected()
    _auto = masktool.source("auto")     # None while masktool.AUTO_MASKS is off
    doc.template_variables.update(
        css=page_css(),
        js_scripts=JS_SCRIPTS, js_version=js_version(),
        main_empty=not workspace.spectra.objects,  # bake the empty-state class into the first paint
        objects_html=workspace.objects_panel_html(),
        processes_html=panel.body_html,
        processes_usage_html=panel.usage_html,
        inspector_html=inspector.info_html(),
        footer_html=workspace.footer_html(),
        swatches_html=swatch_html(sel.color if sel else None),
        edit_html=inspector.edit_html(),
        edit_menu_html=edit_menu.html(),
        style_html=inspector.style_form.html(),
        settings_html=settings.html(),
        synth_run_html=synth_modal.run_form.html(),
        synth_context_html=synth_modal.context_html(),
        synth_actions_html=synth_modal.actions_html(),
        synth_abund_html=synth_modal.abund_form.html(),
        synth_mode_html=synth_modal.mode_form.html(),
        synth_prefs_html=prefs_page.pref_forms["synthprefs"].html(),
        synth_engine_html=prefs_page.pref_forms["synthengine"].html(),
        pref_picker_results=prefs_page.pref_picker_results(),
        pref_labels=PREF_LABELS,
        pref_missing=set(prefs_page.missing),
        pref_show_slurm=load_prefs().get("workflow_runner") == "slurm",
        grid_list_html=prefs_page.grids.list_html(),
        share_list_html={key: share.list_html() for key, share in prefs_page.shares.items()},
        refstars_html=refstars_page.list_html(),
        grid_results_name=prefs_page.grids.results_name,
        synth_grid_html=synth_modal.grid_ui_html(),
        synth_show_file=synth_modal.setup.show_file(),
        synth_file_results_name=synth_modal.atmos_picker.results_name,
        synth_cache_html=prefs_page.cache_html(),
        synth_warning=prefs_page.has_warning,
        synth_atoms_html=synth_atoms.chips_html(),
        synth_is_nlte=synth_modal.is_nlte,
        atom_results_name=synth_atoms.results_name,
        fit_html=synth_modal.fit.html(),
        annot_label=ctx.annot.button_html(),
        linesettings_html=lines.form.html(),
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
        dep1_html=departuretool.form.html(),
        dep_params_html=dep_setup.run_form.html(),
        dep_base_html=dep_setup.base_form.html(),
        dep_grid_html=dep_setup.grid_ui_html(),
        dep_show_file=dep_setup.show_file(),
        dep_file_results_name=dep_setup.atmos_picker.results_name,
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
        results_name=opener.picker.results_name,
    )
    return ctx        # local harnesses drive the document through its context
