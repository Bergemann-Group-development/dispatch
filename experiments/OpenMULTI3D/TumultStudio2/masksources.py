"""The pluggable ways to make a line-mask set.

A :class:`LineMaskSet` can be born in more than one way — read from a line-mask file, started
empty and drawn on the spectrum, or generated automatically from a reference spectrum plus a
line list. Those are not three features but one: each is a *source* that produces the same
object, which then behaves identically everywhere downstream. This module holds the
protocol and the sources that need no synthesis.

A source owns its own widgets and wiring and calls ``on_created(source, mset, message)`` once it
has produced a set; the owning :class:`masktool.MaskTool` decides what happens next (dismiss the
modal, open the view, notify). Nothing about the protocol is synchronous, so a future source that
has to synthesise first can emit from a background job without changing anything here.

Each source's *pane* is declared in the page template (a picker has to be embedded as a Bokeh
root, which cannot be generated from Python), keyed by :attr:`MaskSource.key`; a source's job is
the behaviour behind that pane, not its markup.
"""
from __future__ import annotations

import html
import logging
from pathlib import Path

from tumult import linediag
from controls import Field, Form
from jobs import ResultJobs
from filepicker import FilePicker
from tumult.linemask import is_element


class MaskSource:
    """One way of producing a :class:`linemask.LineMaskSet`.

    Subclasses set the identity attributes, expose whatever Bokeh models and document roots
    they own, and call :meth:`_created` when a set exists.
    """

    key = ""        # pane id in the New line mask modal ("file", "blank", …)
    label = ""      # rail button text
    icon = ""       # tabler icon name for the rail button

    def __init__(self, store, on_created):
        self.store = store              # the shared MaskStore
        self._on_created = on_created   # callback(source, mset, message)

    # --- doc wiring (overridden by sources that own widgets) ------------------------------
    @property
    def models(self):
        """Hidden Bokeh models to collect into the document's bridge column."""
        return []

    @property
    def roots(self):
        """Layouts this source needs added as document roots (e.g. a picker)."""
        return []

    # --- protocol -------------------------------------------------------------------------
    def submit(self, arg=""):
        """The pane's action button was pressed with ``arg``. Sources that act on their own
        (a picker fires on click) leave this alone."""

    def _created(self, mset, message=""):
        self._on_created(self, mset, message)


class FileSource(MaskSource):
    """Read a set from a line-mask table (``input/linemasks``) via an embedded file picker."""

    key, label, icon = "file", "Import from file", "file-import"

    def __init__(self, store, on_created, *, input_dir):
        super().__init__(store, on_created)
        # the read is synchronous, so the picker itself owns both indicators: the clicked
        # row spins and the footer bar is raised for as long as `_pick` runs
        self.picker = FilePicker(self._pick, container_id="linemask-results",
                                 root=input_dir, start_dir=input_dir,
                                 on_busy=lambda on: store.ctx.busy(on))
        self.picker.layout.name = "linemask_picker"

    @property
    def roots(self):
        return [self.picker.layout]

    @property
    def results_name(self):
        """The picker's results-bridge name, so the modal can pre-render on open."""
        return self.picker.results_name

    def _pick(self, path):
        try:
            mset = self.store.load_file(path)
        except Exception as exc:
            logging.exception("Line-mask import failed")
            self._on_created(self, None, f"Could not read line masks: {exc}")
            return
        self._created(mset, f"Imported {len(mset)} line masks from {Path(path).name}.")


class BlankSource(MaskSource):
    """Start an empty set; its windows are drawn on the spectrum."""

    key, label, icon = "blank", "Empty set", "square-plus"

    DEFAULT_NAME = "New line mask"

    def submit(self, arg=""):
        name = str(arg).strip() or self.DEFAULT_NAME
        mset = self.store.add_set(name)
        self._created(mset, f"Started an empty line mask “{mset.name}”.")


class AutoSource(MaskSource):
    """Generate a set by running the diagnostic line selection over an observation.

    Currently unwired: :data:`masktool.AUTO_MASKS` is off, so this source is not registered
    and its pane is not rendered. The class is kept intact — nothing about it is wrong, it is
    the *view* that has not been reconciled: its one-window-per-line output wants overlapping
    windows merged before it reads like an imported set.

    Takes a loaded spectrum and a line list, pulls every transition of the chosen element
    inside the spectrum's range, and derives a window for each
    (:func:`linediag.generate_masks`), measuring only to decide which windows survive — no
    synthesis and no atmosphere, so it needs nothing but the spectrum and the list.

    The run is slow enough to want a worker thread, so ``submit`` hands it to the context's
    background submitter and the set arrives through the same ``on_created`` landing as the
    instant sources; nothing about the protocol had to change to accommodate it.
    """

    key, label, icon = "auto", "Auto-generate", "wand"

    DEFAULTS = {"element": "Fe", "profile": "gauss", "ew_min": 5.0,
                "snr": 100.0, "cont_pct": 0.5}

    def __init__(self, store, on_created, *, ctx, list_observations, get_observation,
                 master_path, profiles):
        super().__init__(store, on_created)
        self.ctx = ctx
        self.list_observations = list_observations   # () -> [(uid, name)]
        self.get_observation = get_observation       # uid -> Spectrum | None
        self.linelist_path = str(master_path)        # defaults to the master line list pref
        self.ref_uid = None
        self.jobs = ResultJobs(ctx)                  # settles this run's own outcome

        self.form = Form("linemaskauto", [
            Field("element", "Element", kind="text", value=self.DEFAULTS["element"],
                  placeholder="e.g. Fe"),
            Field("profile", "EW profile", kind="select", value=self.DEFAULTS["profile"],
                  options=list(profiles)),
            Field("ew_min", "Min EW [mÅ]", kind="number", value=self.DEFAULTS["ew_min"], min=0),
            Field("snr", "SNR (observation)", kind="number", value=self.DEFAULTS["snr"], min=0),
            Field("cont_pct", "Continuum uncertainty [%]", kind="number",
                  value=self.DEFAULTS["cont_pct"], min=0, step=0.1),
        ])
        # the line list is picked with its own inline picker (it is a path, not a form value)
        self.picker = FilePicker(self._set_linelist, container_id="linemask-auto-results",
                                 start_dir=Path(self.linelist_path).parent)
        self.picker.layout.name = "linemask_auto_picker"
        self.path_bridge = ctx.bridges.html("linemask-auto-path")
        self.ref_bridge = ctx.bridges.html("linemask-auto-ref")
        self.open_bridge = ctx.bridges.cls("#linemask-auto-picker", "open")

    # --- doc wiring -----------------------------------------------------------------------
    @property
    def models(self):
        return [*self.form.models]

    @property
    def roots(self):
        return [self.picker.layout]

    @property
    def results_name(self):
        return self.picker.results_name

    # --- panes ----------------------------------------------------------------------------
    def refresh(self):
        """Repaint the observation selector and the line-list path (the canvas changed)."""
        self.ref_bridge.value = self.ref_html()
        self.path_bridge.value = html.escape(Path(self.linelist_path).name)

    def ref_html(self):
        obs = self.list_observations()
        cur = self.ref_uid or (obs[0][0] if obs else None)
        opts = "".join(
            f'<option value="{html.escape(str(uid))}"'
            f'{" selected" if uid == cur else ""}>{html.escape(name)}</option>'
            for uid, name in obs) or '<option value="">— no spectrum loaded —</option>'
        return f'<select class="field-input lm-auto-ref">{opts}</select>'

    def set_ref(self, uid):
        self.ref_uid = uid or None

    def _set_linelist(self, path):
        from htmlbridge import pulse
        self.linelist_path = str(path)
        self.path_bridge.value = html.escape(Path(path).name)
        pulse(self.open_bridge, "off")      # collapse the inline picker

    # --- the run --------------------------------------------------------------------------
    def submit(self, arg=""):
        obs = self.list_observations()
        ref_uid = self.ref_uid or (obs[0][0] if obs else None)
        ref = self.get_observation(ref_uid) if ref_uid else None
        if ref is None:
            self._on_created(self, None, "Load a spectrum to run the diagnostic selection on.")
            return
        v = self.form.values()
        element = str(v["element"] or "").strip().capitalize()
        if not is_element(element):
            self._on_created(self, None, f"“{v['element']}” is not an element symbol.")
            return
        if not Path(self.linelist_path).is_file():
            self._on_created(self, None, f"No line list at {self.linelist_path}.")
            return
        name = dict(obs).get(ref_uid, "spectrum")
        settings = dict(profile=str(v["profile"]), ew_min=float(v["ew_min"]),
                        snr=float(v["snr"]) or None, cont_pct=float(v["cont_pct"]))
        self.jobs.submit(f"linemask-auto:{element}",
                         lambda _working, _token, _progress: self._run(
                             ref, element, name, settings),
                         self._done, on_error=self._failed)

    def _run(self, ref, element, ref_name, settings):
        """Worker thread: candidates of ``element`` in the reference's range, measured."""
        return linediag.generate_mask_set(
            ref, element, linelist=self.linelist_path, reference_name=ref_name,
            pad=self.store.pad, **settings)

    def _failed(self, exc):
        logging.exception("Diagnostic line selection failed", exc_info=exc)
        self._on_created(self, None, str(exc))

    def _done(self, mset):
        if mset is None:
            self._on_created(self, None, "Diagnostic selection found nothing.")
            return
        kept = sum(1 for m in mset.masks if m.include)
        self.store.adopt(mset)
        self._created(mset, f"{mset.element}: {kept} of {len(mset.masks)} candidate windows "
                            f"kept — review the flags and toggles.")
