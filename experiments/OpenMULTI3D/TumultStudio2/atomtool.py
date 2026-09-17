"""The model-atom objects: **Open → Model atom** loads a MULTI-format atom into the
left sidebar; focusing it opens the atom view, two tables in a tab arrangement (the energy
levels and the bound-bound transitions) in place of the spectrum canvas.

The atom is the entry point for NLTE abundance corrections: its right-sidebar panel holds the
**NLTE corrections…** action, which opens the correction modal (:mod:`abundcorrtool`) with
the atom fixed. The correction results computed from an atom are listed indented under its
row, and removing the atom removes them.

The atom is read once with :class:`tumult.atom.AtomClass` (levels, bound-bound and
bound-free records; no collisions, which the tables do not show). Both tables are searched by
a substring query over their text columns and sorted by their headers through the shared
:mod:`tableui` scaffolding.
"""
from __future__ import annotations

import html
import re
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from htmlbridge import parse_action, pulse, toggle
from tableui import SortState, sort_headers
from tool import Tool

TABS = ("levels", "trans")

_ROMAN = ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
# the level records read as "'Level 41 =      z7F2*'": the term is what the tables show
_LEVEL_PREFIX = re.compile(r"^Level\s*\d+\s*=\s*")
# rows rendered at most per table: an Fe atom has >10⁴ lines, and the whole table travels
# through one bridge on every sort or search, so the rest is reached by refining the query
MAX_ROWS = 2000


def term(level):
    """A level's term label without the file's ``Level N =`` prefix and quotes."""
    return _LEVEL_PREFIX.sub("", str(level.label).strip("'\" ")).strip()


def roman(stage):
    """The spectroscopic stage numeral for the atom file's 1-based ion index."""
    stage = int(stage)
    return _ROMAN[stage - 1] if 1 <= stage <= len(_ROMAN) else str(stage)


def _loggf(tr):
    """log gf from the file's oscillator strength and the lower level's weight."""
    gf = float(tr.fosc) * float(tr.lower.g)
    return float(np.log10(gf)) if gf > 0 else float("nan")


@dataclass
class AtomObject:
    """One loaded model atom: the parsed file and its sidebar identity."""
    path: str
    atom: object                     # tumult.atom.AtomClass
    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])

    @property
    def name(self):
        return Path(self.path).name

    @property
    def symbol(self):
        return str(self.atom.id)

    @property
    def lines(self):
        """The bound-bound transitions (the first ``nline`` records of the file)."""
        return self.atom.trans[:self.atom.nline]

    def counts(self):
        a = self.atom
        return int(a.nk), int(a.nline), int(a.ncont)


class AtomTool(Tool):
    """The tool behind the atom objects: the list, the two-table view and the panel."""

    name = "atom"
    kind = "atom"                    # ctx.focus_kind while an atom owns the centre view

    LEVEL_SORT = {"i": True, "label": True, "ion": True, "ecm": True, "g": True, "nbb": False}
    TRANS_SORT = {"kr": True, "lam0": True, "lower": True, "upper": True, "loggf": False,
                  "evi": True, "nq": False}

    def __init__(self, ctx, *, corr):
        Tool.__init__(self, ctx)
        self.corr = corr                 # AbundCorrTool: nested results + the modal
        self.atoms = []                  # list[AtomObject]
        self.selected_uid = None
        self.tab = "levels"
        self.sort = {"levels": SortState(self.LEVEL_SORT, "i"),
                     "trans": SortState(self.TRANS_SORT, "lam0")}
        self.query = {"levels": "", "trans": ""}

        self.own_mode_bridges((".sidebar-right", "atom-active"))
        self.own_panel("atom-panel-content")
        self._view_bridge = ctx.bridges.cls("#main", "atom-mode")
        self._tab_bridges = {t: ctx.bridges.cls("#main", f"atom-tab-{t}") for t in TABS}
        self._objects_bridge = ctx.bridges.html("atom-objects")
        self._tables = {t: ctx.bridges.html(f"atom-{t}-table") for t in TABS}
        self._modal = ctx.bridges.dialog("atom-modal")   # the Open Model Atom picker
        # select / remove / tab / sort / search, space-delimited so a query keeps its spaces
        self.actions = ctx.bridges.action("atomobj-actions")
        self.actions.on_change("value", self._on_action)
        # the wizards hand the single sidebar highlight around through this event
        ctx.subscribe("focus-changed", lambda **_: self.render_objects())

    # --- objects ------------------------------------------------------------------------
    @property
    def models(self):
        return []                        # every bridge is owned through ctx.bridges

    def find(self, uid):
        return next((a for a in self.atoms if a.uid == uid), None)

    def by_path(self, path):
        return next((a for a in self.atoms if a.path == str(path)), None)

    def selected(self):
        return self.find(self.selected_uid)

    def add(self, path, atom, uid=None, select=True):
        """Register a parsed atom; a second load of the same file replaces the first."""
        obj = AtomObject(path=str(path), atom=atom)
        if uid:
            obj.uid = uid
        old = self.by_path(path)
        if old is not None:
            obj.uid = old.uid
            self.atoms[self.atoms.index(old)] = obj
        else:
            self.atoms.append(obj)
        if select:
            pulse(self._modal, "close")      # the picker has done its job
            self.select(obj.uid)
        else:
            self.render_objects()
        return obj

    def select(self, uid):
        obj = self.find(uid)
        if obj is None:
            return
        self.selected_uid = uid
        if self.active:                      # another atom while in the view: repaint
            self.push_panel()
            self.render_tables()
            self.render_objects()
        else:
            self.ctx.tools.activate(self.name)   # -> activate(): the view and the panel

    def remove(self, uid):
        """Drop atom ``uid`` and the correction results computed from it. Returns what
        :meth:`insert` needs to put them back, or ``None``."""
        obj = self.find(uid)
        if obj is None:
            return None
        placed = {"atom": obj, "index": self.atoms.index(obj),
                  "selected": self.selected_uid == uid, "results": []}
        for res in list(self.corr.results_for(obj.path)):
            placed["results"].append(self.corr.remove_result(res.uid))
        self.atoms = [a for a in self.atoms if a.uid != uid]
        if self.selected_uid == uid:
            self.selected_uid = None
            if self.active:
                self.ctx.tools.deactivate_active()
        self.render_objects()
        return placed

    def insert(self, atom, index, selected, results):
        """Put back an atom :meth:`remove` took out, with its correction results at their
        positions; ``selected`` opens its view again."""
        if self.find(atom.uid) is not None:
            raise ValueError(f"Duplicate atom identity: {atom.uid}")
        self.atoms.insert(min(index, len(self.atoms)), atom)
        for placed in reversed(results):      # removed front to back, so back to front
            self.corr.insert_result(**placed)
        if selected:
            self.select(atom.uid)
        self.render_objects()

    def is_focused(self):
        return getattr(self.ctx, "focus_kind", "spectrum") == self.kind

    # --- Tool protocol ------------------------------------------------------------------
    def activate(self):
        obj = self.selected()
        if obj is None:
            return                       # nothing to show: leave ``active`` False
        self.active = True
        self.push_mode("on")
        self.push_panel()
        self.set_tab(self.tab)
        self.render_tables()
        pulse(self._view_bridge, "on")
        self.ctx.focus_kind = self.kind
        self.ctx.emit("focus-changed", kind=self.kind)

    def deactivate(self):
        was = self.active
        self.active = False
        if was:
            pulse(self._view_bridge, "off")
            self.push_mode("off")
            if self.ctx.focus_kind == self.kind:
                self.ctx.focus_kind = "spectrum"
                self.ctx.emit("focus-changed", kind="spectrum")

    def blur(self):
        """Another object took focus: leave the view, keep the selection."""
        if self.active:
            self.ctx.tools.deactivate_active()

    def update(self):
        pass                             # the tables do not depend on the spectra

    # --- actions ------------------------------------------------------------------------
    def _on_action(self, attr, old, new):
        if not new:
            return
        kind, arg = parse_action(new)
        if kind == "select":
            self.select(arg)
        elif kind == "remove":
            self.ctx.removal.request(self.kind, arg)
        elif kind == "tab":
            self.set_tab(arg)
        elif kind == "sort":
            tab, _, col = arg.partition(" ")
            if tab in TABS and self.sort[tab].toggle(col):
                self.render_table(tab)
        elif kind == "search":
            tab, _, query = arg.partition(" ")
            if tab in TABS:
                self.query[tab] = query.strip()
                self.render_table(tab)
        elif kind == "corr":
            obj = self.selected()
            if obj is not None:
                self.corr.open_for(obj)

    def set_tab(self, tab):
        if tab not in TABS:
            return
        self.tab = tab
        for t, bridge in self._tab_bridges.items():
            toggle(bridge, t == tab)

    # --- rendering ----------------------------------------------------------------------
    def render_objects(self):
        self._objects_bridge.value = self.objects_html()

    def objects_html(self):
        """The left-sidebar atoms, each followed by the correction results computed from
        it (indented). Results whose atom is not loaded are listed after the atoms so a
        restored session never hides them."""
        rows = []
        focused = self.is_focused()
        for a in self.atoms:
            sel = " selected" if (focused and a.uid == self.selected_uid) else ""
            nk, nl, nc = a.counts()
            title = f"{a.name} — {nk} levels, {nl} bound-bound, {nc} bound-free"
            results = self.corr.results_for(a.path)
            rows.append(
                f'<div class="object-item atom-item{sel}" data-atom="{a.uid}" '
                f'data-atom-results="{len(results)}" title="{html.escape(title)}">'
                f'<span class="result-badge atom-badge">{html.escape(a.symbol)}</span>'
                f'<span class="object-name">{html.escape(a.name)}</span>'
                f'<button class="result-btn" data-atom-remove="{a.uid}" '
                f'title="Remove from workspace (its correction results go with it)"><i class="ti ti-x"></i></button></div>')
            rows.extend(self.corr.result_rows(results, nested=True))
        orphans = self.corr.results_for(None)
        if orphans:
            rows.append('<div class="section-header result-section-head">Corrections</div>')
            rows.extend(self.corr.result_rows(orphans))
        if not self.atoms and not orphans:
            return ""
        return ('<div class="section-header result-section-head">Model atoms</div>'
                + "".join(rows))

    def panel_html(self):
        obj = self.selected()
        if obj is None:
            return '<div class="hint">Select a model atom.</div>'
        a = obj.atom
        nk, nl, nc = obj.counts()
        results = self.corr.results_for(obj.path)
        kind = "molecule" if a.molec else "atom"
        rows = [
            f'<div class="lm-title"><span class="section-header">{html.escape(obj.symbol)} '
            f'model {kind}</span></div>',
            f'<div class="hint atom-path" title="{html.escape(obj.path)}">'
            f'{html.escape(obj.name)}</div>',
            f'<div class="hint">{nk} levels · {nl} bound-bound · {nc} bound-free '
            f'transitions.<br>A({html.escape(obj.symbol)}) = {float(a.abnd):.2f} · '
            f'mass {float(a.ai):.3f} u.</div>',
            '<button class="text-btn atom-corr">NLTE corrections…</button>',
            '<div class="hint">Synthesises this atom\'s lines in 1D LTE on the reference grid, '
            'measures their equivalent widths, and inverts an NLTE curve of growth on the '
            'chosen atmosphere to the abundance reproducing each one: ΔA = A(NLTE) − A(LTE), '
            'line by line.</div>',
        ]
        if results:
            n = len(results)
            rows.append(f'<div class="hint">{n} correction result{"" if n == 1 else "s"} '
                        'listed under this atom in the left sidebar.</div>')
        return "".join(rows)

    def render_tables(self):
        for tab in TABS:
            self.render_table(tab)

    def render_table(self, tab):
        obj = self.selected()
        if obj is None:
            self._tables[tab].value = ""
            return
        self._tables[tab].value = (self.levels_html(obj) if tab == "levels"
                                   else self.trans_html(obj))

    @staticmethod
    def _matches(query, *texts):
        """Every whitespace-separated token of ``query`` occurs in the joined ``texts``."""
        hay = " ".join(str(t) for t in texts).lower()
        return all(tok in hay for tok in query.lower().split())

    def levels_html(self, obj):
        a = obj.atom
        nbb = np.zeros(a.nk, dtype=int)
        for tr in obj.lines:
            nbb[tr.i] += 1
            nbb[tr.j] += 1
        ground = set(int(i) for i in a.ground_states) | {0}
        levels = list(a.level)
        keys = {"i": lambda lv: lv.i, "label": term, "ion": lambda lv: lv.ion,
                "ecm": lambda lv: lv.ecm, "g": lambda lv: lv.g, "nbb": lambda lv: int(nbb[lv.i])}
        state = self.sort["levels"]
        query = self.query["levels"]
        rows = []
        for i in state.order(levels, keys[state.col]):
            lv = levels[i]
            stage = f"{obj.symbol} {roman(lv.ion)}"
            label = term(lv)
            if query and not self._matches(query, label, stage, lv.i + 1):
                continue
            cls = " atom-ground" if lv.i in ground else ""
            rows.append(
                f'<tr class="atom-row{cls}" data-level="{lv.i}">'
                f'<td>{lv.i + 1}</td><td class="atom-label">{html.escape(label)}</td>'
                f'<td>{html.escape(stage)}</td><td>{lv.ecm:.3f}</td><td>{lv.ev:.4f}</td>'
                f'<td>{lv.g:g}</td><td>{int(nbb[lv.i])}</td></tr>')
        head = "<tr>" + sort_headers(state, [
            ("#", "i"), ("Label", "label"), ("Stage", "ion"), ("E [cm⁻¹]", "ecm"),
            ("E [eV]", None), ("g", "g"), ("b-b lines", "nbb")]) + "</tr>"
        return self._table("atom-levels", head, rows, len(levels))

    def trans_html(self, obj):
        lines = obj.lines
        keys = {"kr": lambda tr: tr.kr, "lam0": lambda tr: tr.lam0,
                "lower": lambda tr: term(tr.lower), "upper": lambda tr: term(tr.upper),
                "loggf": _loggf, "evi": lambda tr: tr.evi, "nq": lambda tr: tr.nq}
        state = self.sort["trans"]
        query = self.query["trans"]
        rows = []
        for i in state.order(lines, keys[state.col]):
            tr = lines[i]
            lam = f"{tr.lam0:.3f}"
            lo, up = term(tr.lower), term(tr.upper)
            if query and not self._matches(query, lam, lo, up, tr.kr + 1):
                continue
            wide = "wide" if tr.iw > 0 else ""
            prof = html.escape(str(getattr(tr, "profile", "")))
            rows.append(
                f'<tr class="atom-row" data-kr="{tr.kr}">'
                f'<td>{tr.kr + 1}</td><td>{lam}</td>'
                f'<td class="atom-label">{html.escape(lo)}</td>'
                f'<td class="atom-label">{html.escape(up)}</td>'
                f'<td>{_loggf(tr):.3f}</td><td>{tr.evi:.3f}</td><td>{tr.nq}</td>'
                f'<td>{wide}</td><td>{prof}</td></tr>')
        head = "<tr>" + sort_headers(state, [
            ("#", "kr"), ("λ air [Å]", "lam0"), ("Lower", "lower"), ("Upper", "upper"),
            ("log gf", "loggf"), ("E_low [eV]", "evi"), ("nq", "nq"), ("", None),
            ("Profile", None)]) + "</tr>"
        return self._table("atom-trans", head, rows, len(lines))

    @staticmethod
    def _table(cls, head, rows, total):
        if not rows:
            return '<div class="linelist-empty">No matches</div>'
        shown = rows[:MAX_ROWS]
        if len(rows) > MAX_ROWS:
            count = (f'<div class="atom-count">First {MAX_ROWS} of {len(rows)} — '
                     'search to narrow the table</div>')
        elif len(rows) != total:
            count = f'<div class="atom-count">{len(rows)} of {total}</div>'
        else:
            count = ""
        return (f'<table class="linelist-table atom-table {cls}"><thead>{head}</thead>'
                f'<tbody>{"".join(shown)}</tbody></table>{count}')
