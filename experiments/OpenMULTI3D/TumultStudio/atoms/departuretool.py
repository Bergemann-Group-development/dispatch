"""The departure-coefficient objects: solving one atmosphere's NLTE populations and
looking at them.

The entry point is a model-atom object (:mod:`atomtool`): its panel's **Departure
coefficients…** action opens this tool's modal with the atom fixed. The modal holds the
atmosphere and stellar parameters (a :class:`synthsetup.SynthSetup` cluster) and the
wavelength-quadrature switch; **Solve** submits one durable ``departures.solve`` workflow
(:mod:`tumult.synthesis.tasks`) on that request. Windows and the line list do not enter a
departure solve, so the modal has neither. The atmosphere is one model: its cluster is
built with ``allow_interpolation=False``, so a registered grid is picked from rather than
interpolated between, since interpolating one would be a departure *grid*
(:class:`tumult.departures.DepartureGrid`), one solve per node, which is a different
object. :meth:`DepartureTool.run` refuses a grid that reaches it another way.

The result object is listed under its atom, is path-backed (a session stores the run
folder and reopens it with :meth:`tumult.departures.Departures.load`) and carries the
populations nowhere: they are read from the run on demand
(:meth:`~tumult.departures.Departures.profile`, off the document thread) and cached per
object for as long as it is listed.

Its view replaces the canvas with the usual top/bottom split: the departure coefficients
against log τ₅₀₀ above, the atom's energy levels or bound-bound transitions below (the
two tables of the atom view, :func:`atomtool.levels_table` / :func:`atomtool.trans_table`).
The selected row is what is drawn, and a view opens on the first row of the table as
displayed, so it shows a curve rather than an empty plot:

* a **level** gives b_i and b_i / b_ground, its departure relative to the ground state of
  its own ionisation stage — the second separates overionisation of the whole stage (a
  common factor, flat in the ratio) from level-to-level redistribution;
* a **transition** gives b_l and b_u of its two levels and their ratio b_u/b_l, which is
  the line source function in units of the Planck function.

Each curve carries a short key for the strip over the plot and the full description the
hover shows. **Pin** keeps the current selection's curves drawn while another row is looked
at: those go on renderers of their own so the figure's inset legend can name each pinned
row once, instead of crowding the strip. The panel's overlay checkboxes draw the same
selection from another solve of the same atom (its keys tagged with that solve's
parameters), which is how two atmospheres or two quadratures are compared.

Curves follow the spectrum convention: each solve is drawn in a run colour of its own,
kept for the session, and the quantities of one selection differ by dash (b solid, the
ratios dashed and dotted). Pinned curves are the same, dimmed.
"""
from __future__ import annotations

import html
import logging
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from bokeh.models import (ColumnDataSource, HoverTool, Legend, LegendItem,
                          Span)
from bokeh.plotting import figure

from tumult.synthesis import AtmosphereGrid

from ..ui.controls import Field, Form
from ..ui.htmlbridge import pulse, toggle
from ..ui.style import BG, BORDER, DASHES, DASH_CSS, MUTED, RUN_COLORS, style_dark
from ..ui.tableui import SortState
from ..ui.tool import Tool
from ..workspace import naming, store
from .atomtool import (LEVEL_SORT, TRANS_SORT, level_order, levels_table, roman,
                       term, trans_order, trans_table)

log = logging.getLogger("tumultstudio.departures")

TABS = ("levels", "trans")
STEPS = {"1": "levels", "2": "trans"}
DEFAULTS = {"line_grid": True}
# the curve colours: the selection in the accent, what it is referred to in the second
# colour, and their ratio in the third. Overlays and pins reuse them at lower alpha.
PIN_ALPHA = 0.45


@dataclass
class DepartureObject:
    """One departure solve as the workspace lists it.

    ``departures`` is the solved run once it exists (``None`` while it is being solved, and
    for a restored object whose run folder is gone). ``star``, ``atmosphere`` and
    ``species`` describe the calculation from the moment it is submitted, so a running
    solve is already described in the sidebar.
    """
    atom_path: str
    name: str = ""
    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])
    departures: object = None
    species: str = ""
    star: dict = field(default_factory=dict)
    atmosphere: str = ""
    progress: dict = field(default_factory=dict)
    termination: str = ""

    @property
    def folder(self):
        return str(self.departures.folder) if self.departures is not None else ""

    @property
    def solved(self):
        return self.departures is not None

    def summary(self):
        if self.progress:
            return self.progress.get("stage", "solving")
        if not self.solved:
            return self.termination or "not on disk"
        star = self.star or self.departures.star
        if star.get("teff"):
            return naming.star(star.get("teff"), star.get("logg"), star.get("feh"),
                               star.get("vmic"))
        return Path(self.atmosphere).name or "solved"


class DepartureTool(Tool):
    """The tool behind the departure objects: the modal, the list, the view and the panel."""

    name = "departures"
    kind = "departures"            # the workspace focus_kind while this view owns the centre
    prefix = "dep"                 # the DOM / action-channel namespace

    def __init__(self, ctx, *, setup, load_async, use_departures):
        Tool.__init__(self, ctx)
        self.setup = setup                     # SynthSetup: atmosphere + stellar parameters
        self.load_async = load_async           # read populations off the document thread
        self.use_departures = use_departures   # hand a solved folder to the synthesis modal
        self.draft = None                      # CalculationDraft, installed by app.build
        self.on_results = lambda: None         # repaints the atom list these hang under
        self.objects = []
        self.selected_uid = None
        self.tab = "levels"
        self.sort = {"levels": SortState(LEVEL_SORT, "i"),
                     "trans": SortState(TRANS_SORT, "lam0")}
        self.query = {"levels": "", "trans": ""}
        self.row = {"levels": None, "trans": None}   # the selected level / transition index
        self.pins = []                     # [(uid, tab, index)] kept drawn behind the selection
        self.overlays = set()              # uids of other solves drawn with the selection
        self._profiles = {}                # uid -> DepartureProfile
        self._loading = set()              # uids whose populations are being read
        self._colors = {}                  # uid -> run colour the solve is drawn in
        self._atoms = {}                   # atom path -> AtomClass
        self._live = {}                    # uid -> the job solving it
        self._pending_atom = None          # the AtomObject the open modal solves for

        vals = {**DEFAULTS, **dict(store.get_section("departures"))}
        self.form = Form("dep1", [
            Field("line_grid", "Line grid from opacity profiles", kind="checkbox",
                  value=bool(vals.get("line_grid", True))),
        ])

        self.own_mode_bridges((".sidebar-right", "dep-active"))
        self.own_panel("dep-panel-content")
        self._view_bridge = ctx.bridges.cls("#main", "dep-mode")
        self._tab_bridges = {t: ctx.bridges.cls("#main", f"dep-tab-{t}") for t in TABS}
        self._tables = {t: ctx.bridges.html(f"dep-{t}-table") for t in TABS}
        self._legend = ctx.bridges.html("dep-legend")
        self._modal = ctx.bridges.dialog("dep-modal", name="dep-modal-bridge")
        self._context_bridge = ctx.bridges.html("dep-context")
        self.actions = ctx.bridges.action("dep-actions")
        self.actions.on_change("value", self._on_action)
        ctx.subscribe("focus-changed", lambda **_: self.on_results())
        ctx.task_importers["departures.solve"] = self.import_task_result

        self.fig = figure(sizing_mode="stretch_both", output_backend="webgl",
                          x_axis_label="log τ₅₀₀", y_axis_label="b = n / n_LTE",
                          y_axis_type="log", name="dep_fig")
        style_dark(self.fig)
        self.src = ColumnDataSource(dict(xs=[], ys=[], color=[], dash=[], width=[],
                                         alpha=[], label=[]))
        self.lines = self.fig.multi_line(xs="xs", ys="ys", line_color="color",
                                         line_dash="dash", line_width="width",
                                         line_alpha="alpha", source=self.src)
        # the hover names a curve in full, which is why the legend keys can stay short
        self.hover = HoverTool(renderers=[self.lines], line_policy="nearest",
                               tooltips=[("", "@label")])
        self.fig.add_tools(self.hover)
        unity = Span(location=1.0, dimension="width", line_color=MUTED, line_width=1,
                     line_dash="dotted")
        self.fig.add_layout(unity)
        # pinned selections are drawn on renderers of their own so the figure's own legend
        # can name them inside the plot, one entry per pinned row, instead of crowding the
        # strip above it with the selection's keys
        self._pinned = []                  # [(ColumnDataSource, renderer)], grown as needed
        # top right: departures grow towards the surface, so the left half of the
        # plot is where the curves are
        self.legend = Legend(items=[], location="top_right", visible=False,
                             background_fill_color=BG, background_fill_alpha=0.75,
                             border_line_color=BORDER, label_text_color=MUTED,
                             label_text_font_size="10px", padding=6, spacing=2,
                             margin=6, glyph_height=10, glyph_width=18)
        self.fig.add_layout(self.legend)

    # --- objects -------------------------------------------------------------------------
    @property
    def models(self):
        return [*self.form.models]

    @property
    def figs(self):
        return [self.fig]

    def find(self, uid):
        return next((o for o in self.objects if o.uid == uid), None)

    def selected(self):
        return self.find(self.selected_uid)

    def results_for(self, atom_path):
        """The solves of ``atom_path`` (``None``: those whose atom is not loaded)."""
        if atom_path is None:
            known = set(self._known_atom_paths())
            return [o for o in self.objects if o.atom_path not in known]
        return [o for o in self.objects if o.atom_path == str(atom_path)]

    def _known_atom_paths(self):
        """Installed by app.build: the paths of the loaded model atoms."""
        return []

    def add(self, obj, select=True):
        self.objects.append(obj)
        self.on_results()
        if select:
            self.select(obj.uid)
        return obj

    def select(self, uid):
        obj = self.find(uid)
        if obj is None:
            return
        self.selected_uid = uid
        self.row = {"levels": None, "trans": None}
        if self.active:
            self.push_panel()
            self.render_tables()
            self.render_curves()
            self.on_results()
        else:
            self.ctx.tools.activate(self.name)
        self._ensure_profile(obj)

    def remove(self, uid):
        """Drop the solve ``uid``, cancelling the work still running on it. Returns what
        :meth:`insert` needs to put it back, or ``None``."""
        obj = self.find(uid)
        if obj is None:
            return None
        placed = {"obj": obj, "index": self.objects.index(obj),
                  "selected": self.selected_uid == uid}
        job = self._live.pop(uid, None)
        if job is not None:
            self.ctx.jobs.request_cancel(job.pid)
        self.objects = [o for o in self.objects if o.uid != uid]
        self.pins = [pin for pin in self.pins if pin[0] != uid]
        self.overlays.discard(uid)
        self._profiles.pop(uid, None)
        if self.selected_uid == uid:
            self.selected_uid = None
            if self.active:
                self.ctx.tools.deactivate_active()
        self.on_results()
        return placed

    def insert(self, obj, index, selected):
        if self.find(obj.uid) is not None:
            raise ValueError(f"Duplicate departures identity: {obj.uid}")
        self.objects.insert(min(index, len(self.objects)), obj)
        if selected:
            self.select(obj.uid)
        self.on_results()

    def rename(self, uid, name):
        """Edit → Rename… on a solve: the object's name, kept distinct in the document."""
        obj = self.find(uid)
        if obj is None:
            return
        obj.name = self.ctx.workspace.unique_name(name, own=obj.name)
        self.on_results()
        if self.selected() is obj:
            self.push_panel()

    def is_computing(self, obj):
        return bool(obj.progress)

    def is_focused(self):
        return self.ctx.workspace.focus_kind == self.kind

    # --- Tool protocol --------------------------------------------------------------------
    def activate(self):
        obj = self.selected()
        if obj is None:
            return
        self.active = True
        self.push_mode("on")
        self.push_panel()
        self.set_tab(self.tab)
        self.render_tables()
        self.render_curves()
        pulse(self._view_bridge, "on")
        self.ctx.workspace.focus(self.kind)

    def deactivate(self):
        was = self.active
        self.active = False
        if was:
            pulse(self._view_bridge, "off")
            self.push_mode("off")
            if self.ctx.workspace.focus_kind == self.kind:
                self.ctx.workspace.focus("spectrum")

    def blur(self):
        if self.active:
            self.ctx.tools.deactivate_active()

    def update(self):
        pass                 # the populations do not depend on the spectra

    def set_tab(self, tab):
        if tab not in TABS:
            return
        self.tab = tab
        for name, bridge in self._tab_bridges.items():
            toggle(bridge, name == tab)

    # --- actions --------------------------------------------------------------------------
    def _on_action(self, attr, old, new):
        # "<kind>:<arg>:<nonce>", the analysis binder's convention; a search query may
        # itself contain a colon, so only the first and last parts are fixed
        if not new:
            return
        parts = new.split(":")
        kind = parts[0]
        arg = ":".join(parts[1:-1]) if len(parts) > 2 else (parts[1] if len(parts) > 1 else "")
        if kind == "select":
            self.select(arg)
        elif kind == "remove":
            self.ctx.removal.request(self.kind, arg)
        elif kind == "step":
            self.set_tab(STEPS.get(arg, "levels"))
            self.render_curves()
            self.push_panel()
        elif kind == "sort":
            if self.sort[self.tab].toggle(arg):
                self.render_table(self.tab)
        elif kind == "search":
            tab, _, query = arg.partition(" ")
            if tab in TABS:
                self.query[tab] = query.strip()
                self.render_table(tab)      # may move the selection onto a shown row
                if tab == self.tab:
                    self.render_curves()
        elif kind == "jump-line":
            self.pick_row(arg)
        elif kind == "pin":
            self.pin()
        elif kind == "unpin":
            self.pins = []
            self.render_curves()
        elif kind == "overlay":
            self.toggle_overlay(arg)
        elif kind == "use":
            self.hand_to_synthesis()
        elif kind == "run":
            self.run()

    def pick_row(self, index):
        try:
            self.row[self.tab] = int(index)
        except (TypeError, ValueError):
            return
        self.render_curves()

    def pin(self):
        """Keep the current selection's curves drawn while another row is looked at."""
        obj, index = self.selected(), self.row[self.tab]
        if obj is None or index is None:
            return
        pin = (obj.uid, self.tab, index)
        if pin not in self.pins:
            self.pins.append(pin)
        self.render_curves()

    def toggle_overlay(self, uid):
        if uid in self.overlays:
            self.overlays.discard(uid)
        else:
            other = self.find(uid)
            if other is None:
                return
            self.overlays.add(uid)
            self._ensure_profile(other)
        self.push_panel()
        self.render_curves()

    def hand_to_synthesis(self):
        obj = self.selected()
        if obj is None or not obj.solved:
            return
        self.use_departures(obj.folder)

    # --- the modal ------------------------------------------------------------------------
    def open_for(self, atom_obj):
        """Open the modal to solve departures of ``atom_obj`` (the atom panel action)."""
        self._pending_atom = atom_obj
        if self.draft is not None:
            self.draft.begin()
        if not store.get_section(self.setup.section).get("_seeded"):
            self.setup.inherit_from()
            store.update_pref(self.setup.section, "_seeded", True)
        self.setup.render()
        self.form.render()
        self._context_bridge.value = self.context_html()
        pulse(self._modal, "open")

    def context_html(self):
        obj = self._pending_atom
        if obj is None:
            return '<div class="hint">No model atom selected.</div>'
        nk, nl, nc = obj.counts()
        return (f'<div class="hint">Equilibrium departure coefficients of '
                f'<b>{html.escape(obj.symbol)}</b> ({html.escape(obj.name)}) — {nk} levels, '
                f'{nl} bound-bound and {nc} bound-free transitions — on one atmosphere. '
                'The populations are saved with the run and can be reused by any synthesis '
                'of this star, over any window. Windows and the line list do not enter a '
                'departure solve.</div>')

    def run(self):
        """**Solve**: submit the departure solve the modal describes."""
        atom = self._pending_atom
        if atom is None:
            return
        request = self.setup.context()
        if request is None:
            self.notify(self.setup.reason or "Configure the atmosphere.")
            pulse(self._modal, "reject")
            return
        if isinstance(request.atmosphere, AtmosphereGrid):
            self.notify("Choose one atmosphere: pick a model from the grid or an individual "
                        "file. An interpolated grid would solve a departure grid instead.")
            pulse(self._modal, "reject")
            return
        if request.departures is not None:
            request = request.update(departures=None)   # the solve carries its own atom
        line_grid = bool(self.form.value("line_grid"))
        if self.draft is not None:
            self.draft.commit()
        store.update_pref("departures", "line_grid", line_grid)
        pulse(self._modal, "close")
        obj = DepartureObject(
            atom_path=str(atom.path), species=atom.symbol, star=dict(request.star.to_dict()),
            atmosphere=str(getattr(request.atmosphere, "path", "") or ""))
        obj.name = self.ctx.workspace.unique_name(
            naming.departures(atom.symbol, request, self.setup.active_grid()))
        obj.progress = {"stage": f"Solving {atom.symbol} departures", "done": 0, "total": 1,
                        "started": time.time()}
        self.add(obj)
        self.notify(f"Solving {atom.symbol} departure coefficients…")
        # the close patch above only reaches the browser once this callback returns, so the
        # submission — a blocking round trip to the workflow daemon — waits a tick, as the
        # other calculation dialogs do
        kwargs = {"atom": str(atom.path), "line_grid": line_grid}
        self.ctx.schedule(lambda: self.submit(obj, request, kwargs))

    def submit(self, obj, request, kwargs):
        """Queue ``obj``'s solve as a durable workflow and follow it."""
        workspace = self.ctx.generation

        def valid():
            return (not self.ctx.closed and self.ctx.generation is workspace
                    and self.find(obj.uid) is obj)

        def done(result):
            if not valid():
                return
            self._live.pop(obj.uid, None)
            self.settle(obj, result)

        def failed(exc):
            if not valid():
                return
            self._live.pop(obj.uid, None)
            obj.progress, obj.termination = {}, "failed"
            self.on_results()
            if self.selected() is obj:
                self.push_panel()

        def native(fraction):
            if not valid() or not obj.progress:
                return
            if fraction is None:
                obj.progress.pop("run", None)
            else:
                obj.progress["run"] = fraction
            self.on_results()

        # the record carries the object's identity, so a document that recovers it
        # continues that object instead of listing a second one
        job = self.ctx.submit_operation(obj.name, "departures.solve", request, done,
                                        kwargs=kwargs, meta={"uid": obj.uid},
                                        on_error=failed, on_native=native)
        if job is None:
            obj.progress, obj.termination = {}, "failed"
            self.on_results()
            return None
        self._live[obj.uid] = job
        job.token.add_callback(
            lambda: self.ctx.schedule(lambda: self._cancelled(obj)))
        return job

    def _cancelled(self, obj):
        self._live.pop(obj.uid, None)
        if self.find(obj.uid) is not obj or obj.solved:
            return
        obj.progress, obj.termination = {}, "cancelled"
        self.on_results()
        if self.selected() is obj:
            self.push_panel()

    def settle(self, obj, departures):
        """A solve finished: keep its run and show it."""
        obj.departures = departures
        obj.progress, obj.termination = {}, ""
        obj.species = departures.species[0] if departures.species else obj.species
        obj.atmosphere = departures.atmosphere or obj.atmosphere
        obj.star = dict(departures.star or obj.star)
        self._profiles.pop(obj.uid, None)
        self.on_results()
        if self.selected() is obj:
            self.push_panel()
            self.render_tables()
        self._ensure_profile(obj)
        self.notify(f"{obj.species or 'Departure'} coefficients solved: {obj.name}")

    def import_task_result(self, result, handle):
        """A ``departures.solve`` record this document did not submit, or one it submitted
        in an earlier session: the object the record names is continued where it is still
        listed, else the solve is listed as a new one."""
        listed = self.find((getattr(handle, "meta", None) or {}).get("uid"))
        if listed is not None:
            return self.settle(listed, result)
        atoms = result.atom_settings[0] if hasattr(result, "atom_settings") else []
        obj = DepartureObject(atom_path=str(result.atom or (atoms[0] if atoms else "")),
                              species=result.species[0] if result.species else "",
                              star=dict(result.star), atmosphere=result.atmosphere)
        obj.name = self.ctx.workspace.unique_name(handle.label
                                                  or f"{obj.species} departures")
        obj.departures = result
        self.add(obj)

    # --- the populations ------------------------------------------------------------------
    def _ensure_profile(self, obj):
        """Read ``obj``'s populations once, off the document thread, and draw them."""
        if obj is None or not obj.solved or obj.uid in self._profiles or obj.uid in self._loading:
            return
        if not obj.departures.saved:
            return
        self._loading.add(obj.uid)
        departures = obj.departures

        def read(token=None):
            return departures.profile()

        def delivered(profile):
            self._loading.discard(obj.uid)
            if self.find(obj.uid) is not obj:
                return
            self._profiles[obj.uid] = profile
            if self.active:
                self.push_panel()
                self.render_curves()

        future = self.load_async(f"read {obj.name}", read, delivered, path=obj.folder)
        if future is None:
            self._loading.discard(obj.uid)

    def atom_of(self, obj):
        """The parsed model atom of ``obj``: the loaded object where there is one, else the
        file, read once and kept. ``None`` when it cannot be read."""
        path = str(obj.atom_path)
        if path in self._atoms:
            return self._atoms[path]
        loaded = self.atom_object(path)
        if loaded is not None:
            self._atoms[path] = loaded.atom
            return loaded.atom
        try:
            from tumult.atom import AtomClass
            self._atoms[path] = AtomClass(path)
        except Exception as exc:
            log.warning("Cannot read the model atom %s: %s", path, exc)
            self._atoms[path] = None
        return self._atoms[path]

    def atom_object(self, path):
        """Installed by app.build: the loaded :class:`atomtool.AtomObject` of ``path``."""
        return None

    # --- the curves -------------------------------------------------------------------------
    def _ground_of(self, atom, level):
        """The index of the ground state of ``level``'s ionisation stage."""
        stage = int(atom.level[level].ion)
        for index in atom.ground_states:
            if int(atom.level[int(index)].ion) == stage:
                return int(index)
        return 0

    def row_title(self, obj, tab, index):
        """What one selected row is, in a few characters: the level's term or the
        transition's wavelength, tagged with the solve when it is not the selected one."""
        atom = self.atom_of(obj)
        if atom is None or index is None:
            return ""
        if tab == "levels":
            return term(atom.level[index]) if 0 <= index < int(atom.nk) else ""
        transitions = atom.trans[:atom.nline]
        return f"{transitions[index].lam0:.1f} Å" if 0 <= index < len(transitions) else ""

    def tag(self, obj):
        """A solve in a few characters, for telling overlaid and pinned curves apart."""
        star = obj.star or (obj.departures.star if obj.solved else {})
        if star.get("teff") and star.get("logg") is not None:
            return f"T{float(star['teff']):.0f} g{float(star['logg']):.2f}"
        return Path(obj.atmosphere).name or obj.species or obj.name

    def color_of(self, obj):
        """The run colour a solve is drawn in: the first palette colour no other solve
        holds, assigned when the solve is first drawn and kept while the tool lives."""
        color = self._colors.get(obj.uid)
        if color is None:
            taken = set(self._colors.values())
            color = next((c for c in RUN_COLORS if c not in taken),
                         RUN_COLORS[len(self._colors) % len(RUN_COLORS)])
            self._colors[obj.uid] = color
        return color

    def curves(self, obj, tab, index, *, alpha=1.0, tag=""):
        """The curves one selected row stands for, as source rows.

        A level gives b and b/b(ground of its stage); a transition gives the two levels'
        b and their ratio, the line source function in units of the Planck function.

        Each row carries both a short ``key``, which is what the legend strip shows, and
        the full ``label`` the hover names it by: the keys have to stay readable side by
        side, while a hovered curve can say everything about itself.
        """
        profile = self._profiles.get(obj.uid)
        atom = self.atom_of(obj)
        if profile is None or atom is None or index is None:
            return []
        nk = min(profile.nk, int(atom.nk))
        ltau = list(map(float, profile.ltau))
        lead = f"{tag} · " if tag else ""
        color = self.color_of(obj)
        solid, dashed, dotted = DASHES[:3]

        def row(values, key, label, line_dash, width=1.8):
            return dict(xs=ltau, ys=[float(v) for v in values], color=color,
                        dash=line_dash, width=width, alpha=alpha,
                        key=f"{lead}{key}", label=f"{lead}{label}")

        if tab == "levels":
            if not 0 <= index < nk:
                return []
            level = atom.level[index]
            ground = self._ground_of(atom, index)
            stage = f"{obj.species or atom.id} {roman(level.ion)}"
            rows = [row(profile.b[index], "b", f"b — {term(level)}", solid)]
            if ground != index:
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative = profile.b[index] / profile.b[ground]
                rows.append(row(relative, "b / b(gs)", f"b / b(ground {stage})", dashed))
            return rows
        transitions = atom.trans[:atom.nline]
        if not 0 <= index < len(transitions):
            return []
        tr = transitions[index]
        lower, upper = int(tr.i), int(tr.j)
        if max(lower, upper) >= nk:
            return []
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = profile.b[upper] / profile.b[lower]
        return [row(profile.b[lower], "b_l", f"b lower — {term(tr.lower)}", solid),
                row(profile.b[upper], "b_u", f"b upper — {term(tr.upper)}", dashed),
                row(ratio, "b_u / b_l", f"b_u / b_l — {tr.lam0:.3f} Å", dotted)]

    def render_curves(self):
        """Repaint the figure: the selection and the solves overlaid on it, with the pinned
        selections behind them in the plot's own legend."""
        obj = self.selected()
        index = self.row[self.tab]
        rows = []
        for uid in sorted(self.overlays):
            other = self.find(uid)
            if other is not None and other is not obj:
                rows += self.curves(other, self.tab, index, tag=self.tag(other))
        if obj is not None:
            rows += self.curves(obj, self.tab, index)
        self.src.data = {key: [r[key] for r in rows]
                         for key in ("xs", "ys", "color", "dash", "width", "alpha",
                                     "key", "label")}
        self._legend.value = self.legend_html(rows)
        self._render_pins()

    def _render_pins(self):
        """Draw the pinned selections on their own renderers and name them in the figure's
        inset legend: one entry per pinned row, whatever the row draws."""
        used, items = 0, []
        for uid, tab, index in self.pins:
            pinned = self.find(uid)
            if pinned is None:
                continue
            curves = self.curves(pinned, tab, index, alpha=PIN_ALPHA)
            if not curves:
                continue
            renderers = []
            for curve in curves:
                source, renderer = self._pin_renderer(used)
                used += 1
                source.data = dict(x=curve["xs"], y=curve["ys"],
                                   label=[f"pinned · {curve['label']}"] * len(curve["xs"]))
                renderer.glyph.line_color = curve["color"]
                renderer.glyph.line_dash = curve["dash"]
                renderer.visible = True
                renderers.append(renderer)
            title = self.row_title(pinned, tab, index)
            tag = self.tag(pinned)
            label = f"{tag} · {title}" if pinned is not self.selected() else title
            items.append(LegendItem(label=label or tag, renderers=renderers))
        for source, renderer in self._pinned[used:]:
            renderer.visible = False
            source.data = dict(x=[], y=[], label=[])
        self.legend.items = items
        self.legend.visible = bool(items)

    def _pin_renderer(self, index):
        """The ``index``-th pinned-curve renderer, created the first time it is needed.

        Each pinned curve needs a renderer of its own for the legend to name it, and its
        own hover entry so it says what it is like every other curve.
        """
        while len(self._pinned) <= index:
            source = ColumnDataSource(dict(x=[], y=[], label=[]))
            renderer = self.fig.line("x", "y", source=source, line_color=MUTED,
                                     line_width=1.4,
                                     line_alpha=PIN_ALPHA, visible=False)
            self._pinned.append((source, renderer))
            self.hover.renderers = [self.lines] + [r for _s, r in self._pinned]
        return self._pinned[index]

    def legend_html(self, rows):
        if not rows:
            obj = self.selected()
            if obj is not None and obj.uid in self._loading:
                return '<span class="hint">Reading the saved populations…</span>'
            if obj is not None and not obj.solved:
                return '<span class="hint">No saved populations to show.</span>'
            noun = "level" if self.tab == "levels" else "transition"
            return f'<span class="hint">Select a {noun} below to plot its departures.</span>'
        title = ""
        obj = self.selected()
        if obj is not None:
            row = self.row_title(obj, self.tab, self.row[self.tab])
            if row:
                title = f'<span class="dep-row">{html.escape(row)}</span>'
        swatches = []
        for r in rows:
            style = (f"border-bottom:2px {DASH_CSS[r['dash']]} {r['color']};"
                     f"opacity:{r['alpha']:.2f}")
            swatches.append(f'<span class="dep-key"><span class="dep-swatch" style="{style}">'
                            f'</span>{html.escape(str(r["key"]))}</span>')
        return title + "".join(swatches)

    # --- rendering ----------------------------------------------------------------------
    def result_rows(self, objects, nested=False):
        """Sidebar rows for ``objects``; ``nested`` indents them under their atom."""
        focused = self.is_focused()
        nest = " result-nested" if nested else ""
        rows = []
        for obj in objects:
            sel = " selected" if (focused and obj.uid == self.selected_uid) else ""
            summary = obj.summary()
            if obj.progress:
                fraction = float(obj.progress.get("run") or 0.0)
                summary = f"{summary} · {100 * fraction:.0f}%"
            rows.append(
                f'<div class="object-item result-item{nest}{sel}" data-tool="{self.prefix}" '
                f'data-result="{obj.uid}">'
                f'<span class="result-badge dep-badge">b</span>'
                f'<span class="result-meta"><span class="object-name" '
                f'title="{html.escape(obj.name)}">{html.escape(obj.name)}</span>'
                f'<span class="result-summary">{html.escape(summary)}</span></span>'
                f'<button class="result-btn" data-result-remove="{obj.uid}" '
                f'title="Remove from workspace"><i class="ti ti-x"></i></button></div>')
        return rows

    def render_tables(self):
        for tab in TABS:
            self.render_table(tab)

    def render_table(self, tab):
        obj = self.selected()
        atom_obj = None if obj is None else self.atom_object(obj.atom_path)
        if atom_obj is None:
            atom = None if obj is None else self.atom_of(obj)
            atom_obj = _AtomView(obj.atom_path, atom) if atom is not None else None
        if atom_obj is None:
            self._tables[tab].value = ('<div class="linelist-empty">The model atom '
                                       f'{html.escape(Path(obj.atom_path).name) if obj else ""}'
                                       ' cannot be read.</div>')
            return
        order = (level_order if tab == "levels" else trans_order)(
            atom_obj, self.sort[tab], self.query[tab])
        self._ensure_row(tab, order)
        table = levels_table if tab == "levels" else trans_table
        self._tables[tab].value = table(atom_obj, self.sort[tab], self.query[tab],
                                        selected=self.row[tab], row_attr="line")

    def _ensure_row(self, tab, order):
        """Select the first row of the table as displayed when nothing else is selected.

        A view that opened on no row would draw nothing, so the first row of the table —
        in its current sort order and through its current search — is what the curves
        belong to until another row is clicked. A selection the search or a different atom
        has taken off the table is replaced the same way rather than left drawing a row
        that cannot be seen.
        """
        if self.row[tab] is not None and self.row[tab] in order:
            return
        self.row[tab] = order[0] if order else None

    def panel_html(self):
        obj = self.selected()
        if obj is None:
            return '<div class="hint">Select a departure solve.</div>'
        rows = [f'<div class="lm-title"><span class="section-header">'
                f'{html.escape(obj.species or "Departures")}</span></div>',
                f'<div class="hint" title="{html.escape(obj.name)}">'
                f'{html.escape(obj.name)}</div>']
        if obj.progress:
            rows.append('<div class="hint">Solving — follow it in the Processes panel; '
                        'closing the tab leaves it running.</div>')
        elif not obj.solved:
            reason = {"cancelled": "The solve was cancelled.",
                      "failed": "The solve failed; the Processes row says why."}.get(
                          obj.termination, "Its run folder no longer holds saved populations.")
            rows.append(f'<div class="hint hint-warn">{reason}</div>')
        else:
            departures = obj.departures
            star = obj.star or departures.star
            profile = self._profiles.get(obj.uid)
            facts = [f'Atmosphere: {html.escape(Path(departures.atmosphere).name or "—")}']
            if star.get("teff"):
                facts.append("Star: " + html.escape(naming.star(
                    star["teff"], star.get("logg"), star.get("feh"), star.get("vmic"))))
            facts.append("Model atom: "
                         + html.escape(Path(departures.atom).name if departures.atom else "—"))
            facts.append("Quadrature: " + ("prepared line grid" if departures.quadrature
                                           else "the atom's own"))
            rows.append('<div class="hint">' + "<br>".join(facts) + "</div>")
            if profile is not None and profile.columns > 1:
                rows.append(f'<div class="hint">{profile.columns} columns, averaged as '
                            '⟨n⟩/⟨n<sub>LTE</sub>⟩ on the mean depth scale.</div>')
            rows.append(self._convergence_html(obj))
            rows.append(f'<div class="hint folder-path" title="{html.escape(obj.folder)}">'
                        f'{html.escape(obj.folder)}</div>')
            rows.append('<button class="text-btn" data-dep-use="1">Use in synthesis…'
                        '</button>')
            rows.append('<div class="hint">Fills the synthesis modal\'s <b>Saved '
                        'departures</b> with this run: its formal solution then runs on '
                        'these populations, for any window.</div>')
        rows.append(self._overlay_html(obj))
        return "".join(rows)

    def _convergence_html(self, obj):
        """The largest population correction the solve ended on, and where."""
        try:
            levels = obj.departures.convergence_levels()
        except OSError as exc:
            log.warning("Cannot read the convergence report of %s: %s", obj.folder, exc)
            return ""
        if not len(levels):
            return ""
        worst = int(np.argmax(levels))
        return ('<div class="hint">Converged to '
                f'{levels[worst]:.1e} at level {worst + 1} (largest final correction).</div>')

    def _overlay_html(self, obj):
        """The other solves of the same atom, each drawable behind the selection."""
        others = [o for o in self.results_for(obj.atom_path) if o is not obj and o.solved]
        if not others:
            return ""
        rows = ['<div class="section-header">Compare</div>']
        for other in others:
            checked = " checked" if other.uid in self.overlays else ""
            rows.append(f'<label class="dep-overlay"><input type="checkbox" '
                        f'data-dep-overlay="{other.uid}"{checked}>'
                        f'{html.escape(other.name)}</label>')
        rows.append('<div class="hint">Draws the same level or transition from another '
                    'solve of this atom in its own colour.</div>')
        return "".join(rows)

    # --- session persistence ----------------------------------------------------------
    def snapshot(self):
        """The listed solves, as the workspace slot stores them (path-backed)."""
        saved = [{"uid": o.uid, "name": o.name, "atom_path": o.atom_path,
                  "folder": o.folder, "species": o.species}
                 for o in self.objects if o.solved]
        # a solve still running has no run to reopen, so it is not saved and cannot be the
        # slot's selection either
        selected = next((d["uid"] for d in saved if d["uid"] == self.selected_uid), None)
        return {"departures": saved, "selected": selected}

    def restore(self, prepared, selected):
        """Replace the listed solves with a saved workspace's prepared ones (each its saved
        identity and the run it was reopened from)."""
        self.objects = [
            DepartureObject(atom_path=p.snap.atom_path, name=p.snap.name, uid=p.snap.uid,
                            departures=p.departures, species=p.snap.species,
                            star=dict(p.departures.star), atmosphere=p.departures.atmosphere)
            for p in prepared]
        self.selected_uid = selected if any(o.uid == selected for o in self.objects) else None
        self.pins, self.overlays = [], set()
        self._profiles, self._atoms = {}, {}
        self.row = {"levels": None, "trans": None}


@dataclass
class _AtomView:
    """The minimum :func:`atomtool.levels_table` needs of an atom object, for a solve whose
    model atom is not loaded in the workspace."""
    path: str
    atom: object

    @property
    def symbol(self):
        return str(self.atom.id)

    @property
    def lines(self):
        return self.atom.trans[:self.atom.nline]
