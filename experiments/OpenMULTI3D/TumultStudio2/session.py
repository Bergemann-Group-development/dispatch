"""Typed schema for a saved workspace (a session "save slot").

``main.snapshot()`` builds a :class:`Workspace` from the live canvas and stores its
``to_dict()`` in the config file (see ``store.save_session``); ``main.restore_workspace()``
reads it back with :meth:`Workspace.from_dict`. Centralising the shape here — instead of
hand-building and re-reading nested dicts with scattered ``.get()`` defaults — documents
the format in one place and makes a malformed/partial slot degrade to sane defaults rather
than raising. The on-disk representation stays plain JSON (nested dicts/lists).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from typing import Optional


def _kept(cls, d):
    """The items of mapping ``d`` whose keys are declared fields of dataclass ``cls``
    (unknown keys are dropped; missing ones fall back to the field defaults)."""
    names = {f.name for f in fields(cls)}
    return {k: v for k, v in (d or {}).items() if k in names}


@dataclass
class ObjectSnap:
    """One restorable spectrum: where it was loaded from plus its display state."""
    uid: str | None = None
    kind: str = ""            # "file" | "run" (path-backed) | "synth" (data embedded below)
    path: str = ""
    name: str = ""
    # synth runs are interpolated, with no backing file/folder, so their x/y are
    # embedded here (a Spectrum.to_dict()); file/run objects leave this None and reload by path.
    data: Optional[dict] = None
    # Info-tab provenance (synth params); FITS metadata is re-read from the file, not stored.
    meta: Optional[dict] = None
    color: Optional[str] = None
    dash: Optional[str] = None
    width: Optional[float] = None
    alpha: Optional[float] = None
    visible: bool = True
    conv: dict = field(default_factory=dict)         # convolution params (vshift/vdop/…)
    norm_knots: Optional[dict] = None                # {"x": [...], "y": [...]} or None
    already_normalised: bool = False

    @classmethod
    def from_dict(cls, d):
        return cls(**_kept(cls, d))


@dataclass
class FilterSnap:
    """The filter state (which method, its parameters, and whether it is shown)."""
    active: bool = True
    method: str = "brault"
    res: Optional[float] = None      # low-pass cutoff resolution
    fmax: Optional[float] = None     # Brault max frequency
    floor: Optional[float] = None    # Brault noise floor (log power)

    @classmethod
    def from_dict(cls, d):
        return cls(**_kept(cls, d))


@dataclass
class Workspace:
    schema_version: int = 1
    objects: list = field(default_factory=list)   # list[ObjectSnap]
    units: dict = field(default_factory=dict)     # {"x": label, "y": label}
    view: dict = field(default_factory=dict)      # {"main": [lo, hi] | None, "nav": ...}
    selected: int = 0
    filter: FilterSnap = field(default_factory=FilterSnap)
    # loaded model atoms, reloaded from their files: [{"uid": ..., "path": ...}]
    atoms: list = field(default_factory=list)
    # tool results, each an opaque {"results": [...], "selected": uid} built by the tool itself
    # (see the snapshot() of AbundanceTool / AbundCorrTool / IonBalanceTool)
    abundance: dict = field(default_factory=dict)
    abundcorr: dict = field(default_factory=dict)
    ionbal: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, d):
        d = d or {}
        return cls(
            objects=[ObjectSnap.from_dict(o) for o in (d.get("objects") or [])],
            units=dict(d.get("units") or {}),
            view=dict(d.get("view") or {}),
            selected=int(d.get("selected") or 0),
            filter=FilterSnap.from_dict(d.get("filter") or {}),
            atoms=[dict(a) for a in (d.get("atoms") or []) if a.get("path")],
            abundance=dict(d.get("abundance") or {}),
            abundcorr=dict(d.get("abundcorr") or {}),
            ionbal=dict(d.get("ionbal") or {}),
        )

    def to_dict(self):
        return asdict(self)
