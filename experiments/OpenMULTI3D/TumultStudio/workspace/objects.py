"""The live objects of one document: the spectra, their identities and the single focus.

A document holds four kinds of object in its left sidebar: spectra, line mask sets, model
atoms and the analysis results. The spectra are owned here (:class:`Spectra`), the mask sets
by :class:`~..masks.maskstore.MaskStore`, the atoms by :class:`~..atoms.atomtool.AtomTool`
and each result list by its :class:`~..analysis.results.ResultCollection`; the
:class:`WorkspaceObjects` of the document reaches all four, and owns what spans them:

* the single focus (``focus_kind``): which list the highlight rests on, and
  :meth:`focused` / :meth:`focus_object` to read and route it the way a click does;
* unique names across every object and every unfinished job (:meth:`unique_name`);
* the Edit menu's Rename and Duplicate of a spectrum or a mask set;
* the sidebar's object list and the footer, which present the spectra and the focus.

A :class:`SpectrumObject` is data and identity alone: the spectrum as loaded and as shown,
its provenance, processing state and line style. The canvas
(:mod:`..spectrum.canvas`) draws it and keeps the Bokeh renderers; the saved workspace
(:mod:`.session`) records it.
"""
from __future__ import annotations

import html
import uuid
from copy import deepcopy

from ..ui.style import DASHES, DASH_CSS, LINE, RUN_COLORS
from . import naming
from .session import FOCUS_TOOL, RESULT_TOOLS

# the sidebar tools that own the left-sidebar focus while open: a panel beside the canvas
FOCUS_TOOLS = (*RESULT_TOOLS, "linemask", "atom", "departures")


class SpectrumObject:
    """One loaded spectrum: its data, provenance, processing state and style."""

    def __init__(self, spectrum, name, dash, uid, color=LINE):
        self.spectrum = spectrum   # the displayed, continuum-normalised spectrum (~0–1)
        self.raw = spectrum        # the originally loaded spectrum, in its native units
        self.norm_knots = None     # (kx, ky) continuum-spline knots in native units, or None
        self.continuum = None      # last continuum evaluated over raw.xdata (for the overlay)
        self.already_normalised = False  # treat the spectrum as pre-normalised (continuum ≡ 1)
        self.name = name
        self.uid = uid
        self.visible = True
        self.source = None  # (kind, path) where it was loaded from, for session save/restore
        self.meta = None  # provenance for the Info tab: FITS header cards or synth params
        self.filtered = None  # filtered Spectrum (reversible); used when canvas.show_filtered
        self.conv = {}  # convolution params applied to the main view (vshift/vdop/…)
        # the line style, as the canvas draws it and the saved workspace records it
        self.color = color
        self.dash = dash
        self.width = 1
        self.alpha = 1.0


class Spectra:
    """The document's spectrum objects, in sidebar order, and the uid of the selected one."""

    def __init__(self):
        self.objects: list[SpectrumObject] = []
        self.selected_uid = None

    def __iter__(self):
        return iter(self.objects)

    def __len__(self):
        return len(self.objects)

    def find(self, uid):
        return next((o for o in self.objects if o.uid == uid), None)

    def selected(self):
        return self.find(self.selected_uid)

    def listed(self):
        """``(uid, name)`` of every object, for the target and reference selectors."""
        return [(o.uid, o.name) for o in self.objects]

    def new(self, spectrum, name=None, color=None, dash=None, uid=None):
        """A :class:`SpectrumObject` ready to be listed, with the defaults an added object
        gets: the next numbered name, the gray line colour, the first line style not in
        use, a fresh identity. A uid already listed is a bug."""
        name = name or f"Object {len(self.objects) + 1}"
        if dash is None:
            used = {o.dash for o in self.objects}
            dash = next((d for d in DASHES if d not in used),
                        DASHES[len(self.objects) % len(DASHES)])
        if uid is None:
            uid = uuid.uuid4().hex
        if self.find(uid) is not None:
            raise ValueError(f"Duplicate observation identity: {uid}")
        return SpectrumObject(spectrum, name, dash, uid, color=color or LINE)

    def next_run_style(self):
        """Pick a (color, dash) for a tumult run: the first unused non-gray palette
        colour, then iterate line styles once every colour is taken."""
        used = {(o.color, o.dash) for o in self.objects}
        used_colors = {c for c, _ in used}
        for color in RUN_COLORS:
            if color not in used_colors:
                return color, DASHES[0]
        for dash in DASHES:
            for color in RUN_COLORS:
                if (color, dash) not in used:
                    return color, dash
        n = len(self.objects)
        return RUN_COLORS[n % len(RUN_COLORS)], DASHES[n % len(DASHES)]

    def append(self, obj, select):
        """List ``obj``; it is selected when asked, or when it is the only object."""
        first = not self.objects
        self.objects.append(obj)
        if select or first:
            self.selected_uid = obj.uid

    def insert(self, obj, index, select):
        """Put ``obj`` back at ``index`` (the removal Undo); a listed uid is a bug."""
        if self.find(obj.uid) is not None:
            raise ValueError(f"Duplicate observation identity: {obj.uid}")
        first = not self.objects
        self.objects.insert(min(index, len(self.objects)), obj)
        if select or first:
            self.selected_uid = obj.uid

    def remove(self, obj):
        """Drop ``obj``; when it was selected, the first remaining object takes its place."""
        self.objects.remove(obj)
        if self.selected_uid == obj.uid:
            self.selected_uid = self.objects[0].uid if self.objects else None

    def select(self, uid):
        """Select ``uid``; whether it is listed."""
        if self.find(uid) is None:
            return False
        self.selected_uid = uid
        return True

    def rename(self, uid, name):
        """Give an object a new display name; whether it changed."""
        obj = self.find(uid)
        name = str(name).strip()
        if obj is None or not name or name == obj.name:
            return False
        obj.name = name
        return True


# shown in the Objects panel while it holds nothing (and nothing is loading): how to
# get a first spectrum onto the canvas
OBJECTS_EMPTY_HINT = (
    '<div class="object-empty">No spectra yet. Open one with '
    '<b>Open → Spectrum</b> <i class="ti ti-folder-open"></i> above, or '
    '<b>New → Synthesise spectrum</b> to compute one.</div>'
)


class WorkspaceObjects:
    """The live objects and the single focus of one document.

    ``spectra`` exists from the start, so the canvas can be built over it; the other
    collections and views are bound with :meth:`bind` once they exist, before anything
    reads the focus. ``on_change`` is the document's repaint after a rename.
    """

    def __init__(self, ctx):
        self.ctx = ctx
        self.spectra = Spectra()
        # which left-sidebar list owns the single selection: "spectrum", "linemask",
        # "atom" or a result kind (session.FOCUS_KINDS); the views move it with
        # :meth:`focus` and the lists repaint on the "focus-changed" event
        self.focus_kind = "spectrum"
        self.canvas = None       # SpectrumCanvas: draws the spectra
        self.masks = None        # MaskStore: the line mask sets
        self.mask_view = None    # MaskTool: the LineMasks view over them
        self.atoms = None        # AtomTool: the model atoms and their view
        self.departures = None   # DepartureTool: the departure solves and their view
        self.results = None      # tool name -> the result view owning that list
        self.on_change = lambda: None

    def bind(self, canvas, *, masks, mask_view, atoms, results, departures):
        self.canvas = canvas
        self.masks = masks
        self.mask_view = mask_view
        self.atoms = atoms
        self.departures = departures
        self.results = results
        assert tuple(results) == RESULT_TOOLS

    # --- the focus ------------------------------------------------------------------------
    def focus(self, kind):
        """Hand the left-sidebar focus to ``kind`` and tell the lists to repaint."""
        self.focus_kind = kind
        self.ctx.emit("focus-changed", kind=kind)

    def focused(self):
        """``(kind, uid)`` of the object holding the single left-sidebar focus."""
        kind = self.focus_kind
        if kind == "spectrum":
            return kind, self.spectra.selected_uid
        if kind == "linemask":
            return kind, self.masks.current_uid
        if kind == "atom":
            return kind, self.atoms.selected_uid
        if kind == "departures":
            return kind, self.departures.selected_uid
        return kind, self.results[FOCUS_TOOL[kind]].results.selected

    def focus_spectrum(self, uid):
        """A click on a spectrum: leave whichever view owns the sidebar, drop the result
        highlights and select it (the canvas repaints with spectrum focus)."""
        tools = self.ctx.tools
        if tools.active_name() in FOCUS_TOOLS:
            tools.deactivate_active()  # selecting a spectrum leaves the wizard / mask view
        for view in self.results.values():
            view.blur()                # leave the centre view, drop the result highlight
        self.canvas.select_object(uid)

    def focus_object(self, kind, uid):
        """Route a focus to ``(kind, uid)`` the way a click on that object does, so the
        panels it opens are the ones a click opens; ``None`` leaves the plain canvas."""
        if kind == "spectrum":
            if uid is not None:
                self.focus_spectrum(uid)
        elif kind == "linemask":
            self.mask_view.select_set(uid)
        elif kind == "atom":
            self.atoms.select(uid)
        elif kind == "departures":
            if uid is not None:
                self.departures.select(uid)
        elif uid is not None:
            self.results[FOCUS_TOOL[kind]].select_result(uid)

    # --- names ----------------------------------------------------------------------------
    def unique_name(self, name, pending=True, own=None):
        """``name`` made unique among the spectrum and result objects (:func:`naming.unique`);
        ``pending`` also avoids the names of this document's unfinished jobs, so two
        submissions of the same calculation deliver distinguishable objects. ``own`` is the
        current name of an object being renamed, which does not count as taken."""
        taken = {o.name for o in self.spectra}
        taken.update(r.name for w in self.results.values() for r in w.results)
        taken.update(o.name for o in self.departures.objects)
        if pending:
            taken.update(job.name for job in self.ctx.jobs.active)
        taken.discard(own)
        return naming.unique(name, taken)

    # --- the spectra ----------------------------------------------------------------------
    def add_synthetic(self, spec, name, meta=None, conv=None):
        """List a delivered synthesis under a distinct run style. The data is embedded on
        save (source path None); a fit also passes its fitted broadening as ``conv``,
        applied live so the Edit tab shows it and the stored spectrum stays raw."""
        color, dash = self.spectra.next_run_style()  # distinct colour, picked at add time
        obj = self.canvas.add_object(spec, name=self.unique_name(name, pending=False),
                                     color=color, dash=dash, source=("synth", None), meta=meta)
        if conv:
            obj.conv = dict(conv)
            self.canvas._refresh_object(obj)  # the repaint that follows rebinds the Edit tab
        return obj

    def rename_spectrum(self, obj, name):
        if self.spectra.rename(obj.uid, self.unique_name(name, own=obj.name)):
            self.on_change()

    def duplicate_spectrum(self, obj):
        """Edit → Duplicate on a spectrum: ``obj`` again under a unique name, with its own
        copy of the data, provenance, continuum, convolution and style, selected so a
        broadening fit or a continuum edit changes the copy only."""
        canvas = self.canvas
        copy = canvas.add_object(deepcopy(obj.raw), name=self.unique_name(obj.name, pending=False),
                                 color=obj.color, dash=obj.dash, source=obj.source,
                                 meta=deepcopy(obj.meta), select=False)
        copy.norm_knots = deepcopy(obj.norm_knots)
        copy.already_normalised = obj.already_normalised
        copy.filtered = deepcopy(obj.filtered)
        copy.conv = dict(obj.conv)
        canvas.normalise_object(copy)
        canvas._refresh_object(copy)
        canvas.set_style(copy, width=obj.width, alpha=obj.alpha)
        canvas.select_object(copy.uid)
        self.ctx.notify(f"Duplicated as “{copy.name}”.", level="success")

    # --- the mask sets --------------------------------------------------------------------
    def rename_set(self, mset, name):
        self.masks.rename(mset.uid, name)

    def duplicate_set(self, mset):
        copy = self.masks.duplicate(mset.uid)
        if copy is not None:
            self.mask_view.select_set(copy.uid)
            self.ctx.notify(f"Duplicated as “{copy.name}”.", level="success")

    # --- the sidebar list and the footer --------------------------------------------------
    def objects_html(self, active_kind=None):
        """The Objects list: one row per spectrum, the selected one highlighted only while
        the focus is on the spectra (``active_kind`` defaults to the current focus)."""
        active_kind = self.focus_kind if active_kind is None else active_kind
        rows = []
        for obj in self.spectra:
            style = DASH_CSS.get(obj.dash, "dashed")
            name = html.escape(obj.name)
            eye = "eye" if obj.visible else "eye-off"
            dim = " hidden-object" if not obj.visible else ""
            # only one item across the Objects + Abundances lists is highlighted at a time
            sel = " selected" if (obj.uid == self.spectra.selected_uid
                                  and active_kind == "spectrum") else ""
            rows.append(
                f'<div class="object-item{dim}{sel}" data-uid="{obj.uid}">'
                f'<span class="dash-swatch" style="border-top-style:{style};border-top-color:{obj.color}"></span>'
                f'<span class="object-name" title="{name}">{name}</span>'
                f'<button class="object-btn" data-action="toggle" data-uid="{obj.uid}" '
                f'title="Hide / show"><i class="ti ti-{eye}"></i></button>'
                f'<button class="object-btn" data-action="remove" data-uid="{obj.uid}" '
                f'title="Remove from workspace (the file stays on disk)"><i class="ti ti-x"></i></button>'
                f'</div>'
            )
        return "".join(rows)

    def objects_panel_html(self, active_kind=None):
        """The Objects panel body: the loaded objects, or the empty-state hint when there is
        nothing to show. Work still running lives in the Processes panel, so a burst of
        background tasks cannot bury the object list."""
        return self.objects_html(active_kind) or OBJECTS_EMPTY_HINT

    def footer_html(self):
        """The footer's object count and the display range they span."""
        n = len(self.spectra)
        if not n:
            return '<span class="pill">0 objects</span>'
        lo, hi = self.canvas.display_bounds()
        rng = f"{lo:.1f} – {hi:.1f} {self.canvas.xlabel}"
        noun = "object" if n == 1 else "objects"
        return f'<span class="pill">{n} {noun}</span><span class="pill">{rng}</span>'
