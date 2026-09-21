"""The saved-workspace schema (a session "save slot") and its headless preparation.

A slot holds one :class:`Workspace`: every object collection of the GUI (spectra, line-mask
sets, model atoms, the three result lists), the canvas view and filter state, the
per-collection selections and the single sidebar focus. ``app.snapshot()`` builds it from
the live document and stores ``to_dict()`` through ``store.save_session``; a load goes the
other way in three steps, so that nothing touches the live workspace until the whole slot
is known to be usable:

1. :meth:`Workspace.from_dict` validates the version, every field type, the identities and
   the references between sections. Anything wrong raises :class:`SessionError` naming the
   field, and the live workspace is untouched.
2. :func:`prepare` turns the validated snapshot into live objects (spectra read from disk or
   rebuilt from embedded data, parsed atoms, mask sets, result objects) with the workspace
   still untouched. External files that cannot be read are the one deliberate partial-load
   case: those objects are skipped and listed in :attr:`PreparedWorkspace.skipped` with
   their identity, path and reason, the references pointing at them are cleared through
   the explicit ``relink`` hook and recorded in :attr:`PreparedWorkspace.unresolved`,
   and a focus or selection that named a skipped object is cleared rather than moved to a
   different object (``relink``; see :mod:`..analysis.results` for the result sections'
   implementation of it). A malformed embedded object (knots, a result, a mask set) is a
   :class:`SessionError`, since the slot itself is broken.
3. ``persistence.WorkspacePersistence`` commits the prepared workspace as one
   document operation, with the rollback and the load report.

Identities: a spectrum, an atom and a result are identified by their string ``uid``, a mask
set by its integer ``uid`` (the :class:`maskstore.MaskStore` counter). References are by
identity only: a result's ``ref_uid`` names a spectrum, a balance result's ``maskset_uid``
a mask set (as a string), a correction result's ``atom_path`` its atom's file; the focus
names one object of one kind. No name-based matching is performed anywhere.

The on-disk representation is plain JSON. Only :data:`SCHEMA_VERSION` is read; an older
slot is rejected with a clear error rather than migrated.
"""
from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

SCHEMA_VERSION = 2

SPECTRUM_KINDS = ("file", "run", "synth")       # path-backed, path-backed, embedded
# what the single sidebar highlight can rest on: the canvas, the mask list, the atom list,
# a departure solve or one of the three result lists (the workspace focus_kind values)
FOCUS_KINDS = ("spectrum", "linemask", "atom", "departures", "result", "corr", "ionbal")
RESULT_TOOLS = ("abundance", "abundcorr", "ionbal")   # the result-list sections
# focus kind -> the result tool whose list carries it
FOCUS_TOOL = {"result": "abundance", "corr": "abundcorr", "ionbal": "ionbal"}
ATOM_TABS = ("levels", "trans")

log = logging.getLogger("tumultstudio.session")


class SessionError(ValueError):
    """The slot cannot be used: unsupported version, wrong shape, or a broken embedded
    object. Raised before any live state changes; the message names the field."""


# --- validation helpers ------------------------------------------------------------------
def _check(cond, where, what):
    if not cond:
        raise SessionError(f"{where}: {what}")


def _mapping(value, where):
    _check(isinstance(value, dict), where, "expected an object")
    return value


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _typed(d, key, kinds, where, *, default=None, required=False):
    """``d[key]`` when it is an instance of ``kinds`` (``None`` allowed unless required);
    ``default`` when absent. A boolean never passes for a numeric field."""
    if key not in d or d[key] is None:
        _check(not required, where, f"missing {key!r}")
        return default
    value = d[key]
    numeric = kinds in (int, float, (int, float))
    if numeric and isinstance(value, bool):
        raise SessionError(f"{where}.{key}: expected a number, got a boolean")
    _check(isinstance(value, kinds), f"{where}.{key}",
           f"unexpected type {type(value).__name__}")
    return value


def _uid(d, key, where, kind=str):
    value = _typed(d, key, kind, where, required=True)
    _check(kind is not str or value.strip(), where, f"{key!r} is empty")
    _check(kind is not int or value > 0, where, f"{key!r} must be positive")
    return value


def _pair(value, where):
    """``[lo, hi]`` of finite numbers with lo < hi, or ``None``."""
    if value is None:
        return None
    _check(isinstance(value, (list, tuple)) and len(value) == 2, where, "expected [lo, hi]")
    lo, hi = value
    _check(_is_number(lo) and _is_number(hi), where, "expected two numbers")
    _check(np.isfinite(lo) and np.isfinite(hi) and lo < hi, where, f"not an interval: {value}")
    return [float(lo), float(hi)]


def _section(d, key):
    """An optional object-valued section: absent (or null) means empty; any other value
    is validated as it stands."""
    return {} if d.get(key) is None else d[key]


def _unique(values, where, what):
    seen = set()
    for v in values:
        _check(v not in seen, where, f"duplicate {what} {v!r}")
        seen.add(v)


# --- sections ----------------------------------------------------------------------------
@dataclass
class SpectrumSnap:
    """One restorable spectrum: where it comes from plus its display state."""
    uid: str
    kind: str                       # one of SPECTRUM_KINDS
    name: str
    path: str = ""                  # file / run objects reload from here
    data: Optional[dict] = None     # synth objects embed a Spectrum.to_dict()
    meta: Optional[dict] = None     # synth provenance (FITS metadata is re-read)
    color: Optional[str] = None
    dash: Optional[str] = None
    width: Optional[float] = None
    alpha: Optional[float] = None
    visible: bool = True
    conv: dict = field(default_factory=dict)
    norm_knots: Optional[dict] = None    # {"x": [...], "y": [...]} of equal length
    already_normalised: bool = False

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        uid = _uid(d, "uid", where)
        where = f"{where}[{uid}]"
        kind = _typed(d, "kind", str, where, required=True)
        _check(kind in SPECTRUM_KINDS, where, f"unknown kind {kind!r}")
        path = _typed(d, "path", str, where, default="")
        data = _typed(d, "data", dict, where)
        if kind == "synth":
            _check(data is not None, where, "a synth object needs embedded data")
        else:
            _check(bool(path), where, f"a {kind} object needs a path")
        knots = _typed(d, "norm_knots", dict, where)
        if knots is not None:
            xs, ys = knots.get("x"), knots.get("y")
            _check(isinstance(xs, list) and isinstance(ys, list) and len(xs) == len(ys),
                   f"{where}.norm_knots", "x and y must be lists of equal length")
            _check(all(_is_number(v) and np.isfinite(v) for v in xs + ys),
                   f"{where}.norm_knots", "knots must be finite numbers")
            knots = {"x": [float(v) for v in xs], "y": [float(v) for v in ys]}
        return cls(
            uid=uid, kind=kind, name=_typed(d, "name", str, where, default=""),
            path=path, data=data, meta=_typed(d, "meta", dict, where),
            color=_typed(d, "color", str, where), dash=_typed(d, "dash", str, where),
            width=_typed(d, "width", (int, float), where),
            alpha=_typed(d, "alpha", (int, float), where),
            visible=_typed(d, "visible", bool, where, default=True),
            conv=dict(_typed(d, "conv", dict, where, default={})),
            norm_knots=knots,
            already_normalised=_typed(d, "already_normalised", bool, where, default=False))


@dataclass
class MaskSnap:
    """The mask store: every set (as ``LineMaskSet.to_dict()``, baseline included so the
    clean/modified state survives), which one is current and which mask is active."""
    sets: list = field(default_factory=list)
    current: Optional[int] = None
    active: Optional[int] = None

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        sets = _typed(d, "sets", list, where, default=[])
        uids = []
        for i, s in enumerate(sets):
            sw = f"{where}.sets[{i}]"
            _mapping(s, sw)
            uids.append(_uid(s, "uid", sw, int))
            masks = _typed(s, "masks", list, sw, default=[])
            for j, m in enumerate(masks):
                _mapping(m, f"{sw}.masks[{j}]")
        _unique(uids, where, "mask set uid")
        current = _typed(d, "current", int, where)
        _check(current is None or current in uids, where, f"current set {current} not saved")
        active = _typed(d, "active", int, where)
        _check(active is None or active >= 0, where, "active must be a mask index")
        return cls(sets=[dict(s) for s in sets], current=current, active=active)


@dataclass
class AtomSnap:
    uid: str
    path: str

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        uid = _uid(d, "uid", where)
        return cls(uid=uid, path=_uid(d, "path", f"{where}[{uid}]"))


@dataclass
class DepartureSnap:
    """One departure solve: its identity, its name and the run folder it is read back from.
    The populations themselves stay in the run (:class:`tumult.departures.Departures`)."""
    uid: str
    name: str
    atom_path: str
    folder: str
    species: str = ""

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        uid = _uid(d, "uid", where)
        return cls(uid=uid, name=_typed(d, "name", str, where, default=""),
                   atom_path=_typed(d, "atom_path", str, where, default=""),
                   folder=_uid(d, "folder", f"{where}[{uid}]"),
                   species=_typed(d, "species", str, where, default=""))


@dataclass
class ResultsSnap:
    """One tool's result list: the serialised results and the selected uid."""
    results: list = field(default_factory=list)
    selected: Optional[str] = None

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        results = _typed(d, "results", list, where, default=[])
        uids = [_uid(_mapping(r, f"{where}.results[{i}]"), "uid", f"{where}.results[{i}]")
                for i, r in enumerate(results)]
        _unique(uids, where, "result uid")
        selected = _typed(d, "selected", str, where)
        _check(selected is None or selected in uids, where,
               f"selected result {selected!r} not saved")
        return cls(results=[dict(r) for r in results], selected=selected)

    @property
    def uids(self):
        return [r["uid"] for r in self.results]


@dataclass
class ViewSnap:
    """Display units and the explicit main/navigation ranges, in display units."""
    xunit: str = ""
    yunit: str = ""
    main_x: Optional[list] = None
    main_y: Optional[list] = None
    nav_x: Optional[list] = None
    nav_y: Optional[list] = None

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        return cls(xunit=_typed(d, "xunit", str, where, default=""),
                   yunit=_typed(d, "yunit", str, where, default=""),
                   **{k: _pair(d.get(k), f"{where}.{k}")
                      for k in ("main_x", "main_y", "nav_x", "nav_y")})


@dataclass
class FilterSnap:
    """The filter state: which method, its parameters, and whether it is shown."""
    active: bool = True
    method: str = "brault"
    res: Optional[float] = None      # low-pass cutoff resolution
    fmax: Optional[float] = None     # Brault max frequency
    floor: Optional[float] = None    # Brault noise floor (log power)

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        return cls(active=_typed(d, "active", bool, where, default=True),
                   method=_typed(d, "method", str, where, default="brault"),
                   res=_typed(d, "res", (int, float), where),
                   fmax=_typed(d, "fmax", (int, float), where),
                   floor=_typed(d, "floor", (int, float), where))


@dataclass
class FocusSnap:
    """The single left-sidebar highlight: which list owns it and which object."""
    kind: str = "spectrum"
    uid: object = None               # str, or int for a mask set; None for no highlight

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        kind = _typed(d, "kind", str, where, default="spectrum")
        _check(kind in FOCUS_KINDS, where, f"unknown focus kind {kind!r}")
        uid = d.get("uid")
        _check(uid is None or isinstance(uid, str) or (isinstance(uid, int)
               and not isinstance(uid, bool)), where, "uid must be a string, an integer or null")
        return cls(kind=kind, uid=uid)


@dataclass
class ToolsSnap:
    """Tool-view state worth carrying: the active tool, the atom view's tab, each result
    wizard's step and its jumped-to (selected diagnostic) line index."""
    active: Optional[str] = None
    atom_tab: str = "levels"
    steps: dict = field(default_factory=dict)      # tool name -> step (int >= 1)
    jumped: dict = field(default_factory=dict)     # tool name -> line index or None

    @classmethod
    def from_dict(cls, d, where):
        d = _mapping(d, where)
        tab = _typed(d, "atom_tab", str, where, default="levels")
        _check(tab in ATOM_TABS, where, f"unknown atom tab {tab!r}")
        steps = _typed(d, "steps", dict, where, default={})
        jumped = _typed(d, "jumped", dict, where, default={})
        for name, v in steps.items():
            _check(name in RESULT_TOOLS and isinstance(v, int) and not isinstance(v, bool)
                   and v >= 1, f"{where}.steps", f"bad entry {name!r}: {v!r}")
        for name, v in jumped.items():
            _check(name in RESULT_TOOLS and (v is None or isinstance(v, int)
                   and not isinstance(v, bool) and v >= 0),
                   f"{where}.jumped", f"bad entry {name!r}: {v!r}")
        return cls(active=_typed(d, "active", str, where), atom_tab=tab,
                   steps=dict(steps), jumped=dict(jumped))


@dataclass
class Workspace:
    schema_version: int = SCHEMA_VERSION
    spectra: list = field(default_factory=list)       # list[SpectrumSnap]
    spectrum_selected: Optional[str] = None
    masks: MaskSnap = field(default_factory=MaskSnap)
    atoms: list = field(default_factory=list)         # list[AtomSnap]
    atom_selected: Optional[str] = None
    departures: list = field(default_factory=list)    # list[DepartureSnap]
    departure_selected: Optional[str] = None
    results: dict = field(default_factory=dict)       # tool name -> ResultsSnap
    view: ViewSnap = field(default_factory=ViewSnap)
    filter: FilterSnap = field(default_factory=FilterSnap)
    focus: FocusSnap = field(default_factory=FocusSnap)
    tools: ToolsSnap = field(default_factory=ToolsSnap)

    @classmethod
    def from_dict(cls, d):
        """Validate a stored snapshot completely; :class:`SessionError` on the first
        problem, which names the offending field."""
        d = _mapping(d, "workspace")
        version = d.get("schema_version")
        if version != SCHEMA_VERSION:
            raise SessionError(f"Saved workspace version {version!r} is not supported "
                               f"(expected {SCHEMA_VERSION})")
        spectra = [SpectrumSnap.from_dict(s, f"spectra[{i}]")
                   for i, s in enumerate(_typed(d, "spectra", list, "workspace", default=[]))]
        _unique([s.uid for s in spectra], "spectra", "spectrum uid")
        atoms = [AtomSnap.from_dict(a, f"atoms[{i}]")
                 for i, a in enumerate(_typed(d, "atoms", list, "workspace", default=[]))]
        _unique([a.uid for a in atoms], "atoms", "atom uid")
        _unique([a.path for a in atoms], "atoms", "atom path")
        departures = [DepartureSnap.from_dict(v, f"departures[{i}]")
                      for i, v in enumerate(_typed(d, "departures", list, "workspace",
                                                   default=[]))]
        _unique([v.uid for v in departures], "departures", "departures uid")
        results_d = _typed(d, "results", dict, "workspace", default={})
        unknown = set(results_d) - set(RESULT_TOOLS)
        _check(not unknown, "results", f"unknown section {sorted(unknown)}")
        results = {name: ResultsSnap.from_dict(_section(results_d, name), f"results.{name}")
                   for name in RESULT_TOOLS}
        ws = cls(
            spectra=spectra,
            spectrum_selected=_typed(d, "spectrum_selected", str, "workspace"),
            masks=MaskSnap.from_dict(_section(d, "masks"), "masks"),
            atoms=atoms,
            atom_selected=_typed(d, "atom_selected", str, "workspace"),
            departures=departures,
            departure_selected=_typed(d, "departure_selected", str, "workspace"),
            results=results,
            view=ViewSnap.from_dict(_section(d, "view"), "view"),
            filter=FilterSnap.from_dict(_section(d, "filter"), "filter"),
            focus=FocusSnap.from_dict(_section(d, "focus"), "focus"),
            tools=ToolsSnap.from_dict(_section(d, "tools"), "tools"))
        # references between the sections
        _check(ws.spectrum_selected is None or ws.spectrum_selected in ws.spectrum_uids,
               "spectrum_selected", f"{ws.spectrum_selected!r} is not a saved spectrum")
        _check(ws.atom_selected is None or ws.atom_selected in ws.atom_uids,
               "atom_selected", f"{ws.atom_selected!r} is not a saved atom")
        _check(ws.departure_selected is None
               or ws.departure_selected in ws.departure_uids,
               "departure_selected",
               f"{ws.departure_selected!r} is not a saved departure solve")
        _check(ws.focus.uid is None or ws.focus.uid in ws.uids_of(ws.focus.kind),
               "focus", f"{ws.focus.kind} {ws.focus.uid!r} is not a saved object")
        return ws

    # --- identities ------------------------------------------------------------------
    @property
    def spectrum_uids(self):
        return [s.uid for s in self.spectra]

    @property
    def atom_uids(self):
        return [a.uid for a in self.atoms]

    @property
    def departure_uids(self):
        return [d.uid for d in self.departures]

    @property
    def maskset_uids(self):
        return [s["uid"] for s in self.masks.sets]

    def uids_of(self, kind):
        """The identities the focus of ``kind`` may name."""
        if kind == "spectrum":
            return self.spectrum_uids
        if kind == "linemask":
            return self.maskset_uids
        if kind == "atom":
            return self.atom_uids
        if kind == "departures":
            return self.departure_uids
        return self.results[FOCUS_TOOL[kind]].uids

    def to_dict(self):
        return asdict(self)


def summary(snapshot):
    """A one-line description of a stored slot for the Sessions list, tolerant of a slot
    this version cannot load (it says so instead of raising)."""
    if not isinstance(snapshot, dict):
        return "unreadable"
    if snapshot.get("schema_version") != SCHEMA_VERSION:
        return f"unsupported format (version {snapshot.get('schema_version')!r})"

    def length(value):
        return len(value) if isinstance(value, list) else 0

    masks = snapshot.get("masks") if isinstance(snapshot.get("masks"), dict) else {}
    sets = masks.get("sets") if isinstance(masks.get("sets"), list) else []
    n_masks = sum(1 for s in sets if isinstance(s, dict) and (s.get("label") or s.get("masks")))
    results = snapshot.get("results") if isinstance(snapshot.get("results"), dict) else {}
    n_results = sum(length((results.get(t) or {}).get("results"))
                    for t in RESULT_TOOLS if isinstance(results.get(t), dict))
    parts = [plural(length(snapshot.get("spectra")), "spectrum", "spectra")]
    for n, noun in ((n_masks, "line mask"), (length(snapshot.get("atoms")), "atom"),
                    (length(snapshot.get("departures")), "departure solve"),
                    (n_results, "result")):
        if n:
            parts.append(plural(n, noun))
    return " · ".join(parts)


def is_blank(snapshot):
    """Whether a workspace dict holds nothing worth keeping: no spectra, atoms, named or
    filled line-mask sets or tool results. A fresh document is blank; the autosave skips
    it so an accidental new tab does not overwrite the last real session."""
    if not isinstance(snapshot, dict):
        return True
    masks = snapshot.get("masks") if isinstance(snapshot.get("masks"), dict) else {}
    sets = masks.get("sets") if isinstance(masks.get("sets"), list) else []
    results = snapshot.get("results") if isinstance(snapshot.get("results"), dict) else {}
    return not (snapshot.get("spectra") or snapshot.get("atoms")
                or snapshot.get("departures")
                or any(isinstance(s, dict) and (s.get("label") or s.get("masks")) for s in sets)
                or any((results.get(t) or {}).get("results") for t in RESULT_TOOLS
                       if isinstance(results.get(t), dict)))


def plural(n, noun, plural=None):
    return f"{n} {noun if n == 1 else (plural or noun + 's')}"


# --- preparation: live objects, off the live workspace ------------------------------
@dataclass
class Skipped:
    kind: str          # "spectrum" | "atom" | "departures" | "workflow"
    uid: str
    name: str
    path: str
    reason: str

    def __str__(self):
        return f"{self.kind} {self.name or self.uid} ({self.path}): {self.reason}"


@dataclass
class PreparedSpectrum:
    snap: SpectrumSnap
    spectrum: object                 # tumult Spectrum
    meta: Optional[dict]


@dataclass
class PreparedAtom:
    snap: AtomSnap
    atom: object                     # tumult.atom.AtomClass


@dataclass
class PreparedDepartures:
    snap: DepartureSnap
    departures: object               # tumult.departures.Departures


@dataclass
class PreparedWorkspace:
    """Everything a commit needs, already built. Also what a rollback re-applies, so the
    live workspace is captured into one of these before a load replaces it."""
    spectra: list = field(default_factory=list)          # list[PreparedSpectrum]
    spectrum_selected: Optional[str] = None
    mask_sets: list = field(default_factory=list)        # list[LineMaskSet]
    mask_current: Optional[int] = None
    mask_active: Optional[int] = None
    atoms: list = field(default_factory=list)            # list[PreparedAtom]
    atom_selected: Optional[str] = None
    departures: list = field(default_factory=list)       # list[PreparedDepartures]
    departure_selected: Optional[str] = None
    results: dict = field(default_factory=dict)          # tool -> (list[result], selected)
    view: ViewSnap = field(default_factory=ViewSnap)
    filter: FilterSnap = field(default_factory=FilterSnap)
    focus: FocusSnap = field(default_factory=FocusSnap)
    tools: ToolsSnap = field(default_factory=ToolsSnap)
    skipped: list = field(default_factory=list)          # list[Skipped]
    unresolved: list = field(default_factory=list)       # list[str]

    def counts(self):
        """Restored objects per collection; a mask set counts only when it would be
        listed in the sidebar (the store's empty scratch set does not)."""
        return {"spectra": len(self.spectra),
                "masks": sum(1 for s in self.mask_sets if s.label or s.masks),
                "atoms": len(self.atoms), "departures": len(self.departures),
                "results": sum(len(r) for r, _ in self.results.values())}


def prepare(ws, *, load_spectrum, load_atom, load_departures, deserialise, relink):
    """Build the live objects of validated ``ws`` without touching the workspace.

    ``load_spectrum(snap) -> Spectrum``, ``load_atom(path) -> AtomClass`` and
    ``load_departures(folder) -> Departures`` read the external files; a failure there skips
    that object (the partial-load policy) and is recorded with its reason. ``deserialise(name, records) -> list`` builds one result
    section's objects and ``relink(name, results, spectrum_uids, maskset_uids, atom_paths)
    -> list[str]`` clears the references that point at objects this load does not restore;
    both are plain functions of the saved records, so a slot is interpreted without
    constructing any view. Every other failure is a :class:`SessionError`.
    """
    from tumult.linemask import LineMaskSet

    prepared = PreparedWorkspace(view=ws.view, filter=ws.filter, tools=ws.tools)
    for snap in ws.spectra:
        try:
            spec = load_spectrum(snap)
        except Exception as exc:
            if snap.kind == "synth":
                raise SessionError(f"spectra[{snap.uid}]: embedded data is unusable: "
                                   f"{exc}") from exc
            log.warning("Skipping %s %r (%s): %s", snap.kind, snap.name, snap.path, exc)
            prepared.skipped.append(Skipped("spectrum", snap.uid, snap.name, snap.path,
                                            f"{type(exc).__name__}: {exc}"))
            continue
        meta = snap.meta if snap.kind == "synth" else getattr(spec, "meta", None)
        prepared.spectra.append(PreparedSpectrum(snap, spec, meta))
    for i, d in enumerate(ws.masks.sets):
        try:
            mset = LineMaskSet.from_dict(d)
        except Exception as exc:
            raise SessionError(f"masks.sets[{i}]: {exc}") from exc
        prepared.mask_sets.append(mset)
    prepared.mask_current, prepared.mask_active = ws.masks.current, ws.masks.active
    for snap in ws.atoms:
        try:
            atom = load_atom(snap.path)
        except Exception as exc:
            log.warning("Skipping atom %s: %s", snap.path, exc)
            prepared.skipped.append(Skipped("atom", snap.uid, Path(snap.path).name, snap.path,
                                            f"{type(exc).__name__}: {exc}"))
            continue
        prepared.atoms.append(PreparedAtom(snap, atom))
    for snap in ws.departures:
        try:
            departures = load_departures(snap.folder)
        except Exception as exc:
            log.warning("Skipping the departure solve %s: %s", snap.folder, exc)
            prepared.skipped.append(Skipped("departures", snap.uid, snap.name, snap.folder,
                                            f"{type(exc).__name__}: {exc}"))
            continue
        prepared.departures.append(PreparedDepartures(snap, departures))

    spectrum_uids = {p.snap.uid for p in prepared.spectra}
    maskset_uids = {s.uid for s in prepared.mask_sets}
    atom_uids = {p.snap.uid for p in prepared.atoms}
    atom_paths = {p.snap.path for p in prepared.atoms}
    for name in RESULT_TOOLS:
        section = ws.results[name]
        try:
            results = deserialise(name, section.results)
        except SessionError:
            raise
        except Exception as exc:
            raise SessionError(f"results.{name}: {exc}") from exc
        prepared.unresolved += relink(name, results, spectrum_uids, maskset_uids,
                                      atom_paths)
        prepared.results[name] = (results, section.selected)

    # selections and the focus follow the skip: a missing object clears them, it never
    # hands the highlight to a different object
    prepared.spectrum_selected = (ws.spectrum_selected
                                  if ws.spectrum_selected in spectrum_uids else None)
    prepared.atom_selected = ws.atom_selected if ws.atom_selected in atom_uids else None
    departure_uids = {p.snap.uid for p in prepared.departures}
    prepared.departure_selected = (ws.departure_selected
                                   if ws.departure_selected in departure_uids else None)
    focus = FocusSnap(ws.focus.kind, ws.focus.uid)
    if (focus.kind == "spectrum" and focus.uid not in spectrum_uids) or (
            focus.kind == "atom" and focus.uid not in atom_uids) or (
            focus.kind == "departures" and focus.uid not in departure_uids):
        focus = FocusSnap("spectrum", None)
    prepared.focus = focus
    return prepared


@dataclass
class LoadReport:
    """What a committed load did: the counts restored, the objects skipped (with their
    identity, path and reason), the references cleared, the jobs the boundary
    cancelled and what became of the slot's workflow references (``"<name>: <outcome>"``;
    a reference whose record is gone is a ``Skipped`` of kind ``workflow``).
    ``message()`` is the footer line, ``details()`` one line per item."""
    name: str
    counts: dict
    skipped: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)
    cancelled_jobs: list = field(default_factory=list)
    workflows: list = field(default_factory=list)

    @property
    def partial(self):
        return bool(self.skipped or self.unresolved)

    def message(self):
        c = self.counts
        parts = [plural(c["spectra"], "spectrum", "spectra"), plural(c["masks"], "line mask"),
                 plural(c["atoms"], "atom"), plural(c["departures"], "departure solve"),
                 plural(c["results"], "result")]
        text = f"Loaded “{self.name}”: " + ", ".join(parts)
        if self.skipped:
            text += "; skipped " + "; ".join(str(s) for s in self.skipped)
        if self.unresolved:
            n = len(self.unresolved)
            text += f"; {n} reference{'' if n == 1 else 's'} cleared (see the log)"
        if self.cancelled_jobs:
            text += f"; cancelled {', '.join(self.cancelled_jobs)}"
        if self.workflows:
            text += f"; workflows: {', '.join(self.workflows)}"
        return text

    def details(self):
        return ([f"Loading “{self.name}” skipped {s}" for s in self.skipped]
                + [f"Loading “{self.name}”: {u}" for u in self.unresolved])
