"""The LineMask sets for TumultStudio.

One :class:`MaskStore` holds every imported :class:`linemask.LineMaskSet` — the app's
line-mask *objects*, listed in the left sidebar next to the spectra — plus which one is
*current*. The current set is the single source of truth every consumer reads: the canvas
draws it, the LineMasks view edits it, and (later) the abundance / synthesize-fit workflows
adopt it. Mutations broadcast the context's ``"masks-changed"`` event so every subscriber
re-syncs, exactly like the existing ``"changed"`` / ``"filter-state"`` events.

The whole single-set API (``masks`` / ``pad`` / ``add`` / ``set_fit`` / …) delegates to the
current set, so callers that only ever cared about "the masks" are unaffected by the
multi-set layering.

The store is deliberately GUI-free (it only touches the context's event bus) so it can be
unit-tested and serialised into a session snapshot on its own.
"""
from __future__ import annotations

import itertools
from copy import deepcopy
from pathlib import Path

from ..workspace import naming

from tumult.linemask import (LineMask, LineMaskSet, parse_linemask_file,
                      write_linemask_file)


class MaskStore:
    def __init__(self, ctx):
        self.ctx = ctx
        # every imported / drawn set, in import order. The store always keeps at least one
        # (a scratch set that box-select can drop a mask into before anything is imported);
        # `is_listed` decides which of them earn a sidebar row.
        self._uids = itertools.count(1)
        self.sets = [self._new_set()]
        self.current_uid = self.sets[0].uid
        self.active = None      # index of the mask being edited on the canvas, or None

    def _new_set(self, mset=None) -> LineMaskSet:
        mset = mset if mset is not None else LineMaskSet()
        mset.uid = next(self._uids)
        return mset

    # --- set-level queries --------------------------------------------------------------
    @property
    def mset(self) -> LineMaskSet:
        """The set currently being viewed / edited (never ``None``)."""
        # `is None`, not a truth test: LineMaskSet defines __len__, so an *empty* set — which
        # is exactly what "Empty line mask" creates — is falsy and would fall through here.
        found = self.find(self.current_uid)
        return self.sets[0] if found is None else found

    def find(self, uid) -> LineMaskSet | None:
        return next((s for s in self.sets if s.uid == uid), None)

    @staticmethod
    def is_listed(mset) -> bool:
        """Does this set deserve a sidebar row? Imported sets always do; the scratch set
        only once something has been drawn into it."""
        return bool(mset.label or mset.masks)

    @property
    def listed(self) -> list:
        """The sets shown as sidebar objects, in import order."""
        return [s for s in self.sets if self.is_listed(s)]

    # --- mask-level queries (delegate to the current set) --------------------------------
    @property
    def masks(self) -> list:
        return self.mset.masks

    @property
    def pad(self) -> float:
        return self.mset.pad

    def active_mask(self) -> LineMask | None:
        if self.active is None or not (0 <= self.active < len(self.masks)):
            return None
        return self.masks[self.active]

    # --- mutations (each ends by broadcasting) ------------------------------------------
    def _changed(self):
        self.ctx.emit("masks-changed")

    def load_file(self, path):
        """Import ``path`` as a new set, make it current, and return it."""
        mset = self._new_set(parse_linemask_file(path, pad=self.mset.pad))
        mset.mark_clean()       # the file on disk is this set's pristine state
        self.sets.append(mset)
        self.current_uid = mset.uid
        self.active = 0 if mset.masks else None
        self._changed()
        return mset

    def add_set(self, label="", element=""):
        """Start a new, empty set under ``label`` and make it current — the "blank sheet"
        counterpart to :meth:`load_file`, for masks drawn from scratch. It has no baseline,
        so it reads as unsaved rather than modified until it is written out."""
        mset = self._new_set(LineMaskSet(pad=self.mset.pad, label=str(label).strip(),
                                         element=str(element).strip()))
        self.sets.append(mset)
        self.current_uid = mset.uid
        self.active = None
        self._changed()
        return mset

    def adopt(self, mset):
        """Take a set built elsewhere (a generator running off the document thread), give it
        an identity and make it current."""
        self.sets.append(self._new_set(mset))
        self.current_uid = mset.uid
        self.active = 0 if mset.masks else None
        self._changed()
        return mset

    def duplicate(self, uid):
        """Add a copy of set ``uid`` under a unique name and make it current. The copy has
        the same windows, element and padding but neither a file nor a baseline, so it is
        unsaved until written out and editing it leaves the original as it is."""
        mset = self.find(uid)
        if mset is None:
            return None
        copy = deepcopy(mset)
        copy.source_path, copy.baseline = "", {}
        copy.label = naming.unique(mset.name, {s.name for s in self.sets})
        return self.adopt(copy)

    def rename(self, uid, name):
        """Give set ``uid`` the label ``name``, made unique among the sets. The label is not
        part of the set's content, so a rename does not mark it modified."""
        mset = self.find(uid)
        if mset is None:
            return None
        mset.label = naming.unique(str(name).strip(),
                                   {s.name for s in self.sets if s is not mset})
        self._changed()
        return mset

    def set_include(self, indices, on):
        """Switch the masks at ``indices`` on/off. The toggle curates what a save writes out
        and what the wizards measure — a flagged candidate stays visible but is left out."""
        indices = [i for i in indices if 0 <= i < len(self.masks)]
        for i in indices:
            self.masks[i].include = bool(on)
        if indices:
            self._changed()

    def set_current(self, uid):
        """Switch which set is viewed / edited (no-op for an unknown uid)."""
        mset = self.find(uid)
        if mset is None:
            return False
        self.current_uid = uid
        self.active = 0 if mset.masks else None
        self._changed()
        return True

    def remove_set(self, uid):
        """Drop a whole set. The store keeps at least one (a fresh scratch set).
        Returns what :meth:`insert_set` needs to put it back, or ``None``."""
        mset = self.find(uid)
        if mset is None:
            return None
        placed = {"mset": mset, "index": self.sets.index(mset),
                  "current": self.current_uid == uid, "active": self.active, "stand_in": None}
        self.sets.remove(mset)
        if not self.sets:
            placed["stand_in"] = self._new_set()
            self.sets.append(placed["stand_in"])
        if self.current_uid == uid:
            self.current_uid = self.sets[-1].uid
            self.active = None
        self._changed()
        return placed

    def insert_set(self, mset, index, current, active, stand_in=None):
        """Put back a set :meth:`remove_set` took out (its returned dict is the arguments),
        with its uid, position and edits; ``current``/``active`` restore it as the viewed
        set. The scratch set that stood in while the store was otherwise empty goes again
        unless something was drawn into it."""
        if self.find(mset.uid) is not None:
            raise ValueError(f"Duplicate line mask set identity: {mset.uid}")
        if stand_in is not None and not self.is_listed(stand_in) and stand_in in self.sets:
            self.sets.remove(stand_in)
        self.sets.insert(min(index, len(self.sets)), mset)
        if current or self.find(self.current_uid) is None:
            self.current_uid, self.active = mset.uid, active
        self._changed()
        return mset

    def revert(self):
        """Restore the current set to the file it was imported from / last saved to."""
        if not self.mset.revert():
            return False
        self.active = 0 if self.masks else None
        self._changed()
        return True

    def save_as(self, path):
        """Write the current set to ``path`` and adopt it as the set's new pristine origin
        (so the edited copy stops reading as modified, and reverts go back to it)."""
        mset = self.mset
        path = Path(path)
        write_linemask_file(path, mset)
        mset.source_path = str(path)
        mset.label = path.stem
        mset.mark_clean()
        self._changed()
        return path

    def default_save_path(self):
        """Where "Save as…" starts: the source file's ``<stem>_edited.txt`` sibling, or an
        export in the working directory for a set drawn from scratch."""
        src = self.mset.source_path
        if src:
            return Path(src).with_name(Path(src).stem + "_edited.txt")
        return Path.cwd() / f"{self.mset.name.replace(' ', '_')}.txt"

    def add(self, lmin, lmax, *, lam0=float("nan"), note="") -> int:
        """Append a new mask over ``[lmin, lmax]`` (synth window derived from ``pad``); select it."""
        mset = self.mset
        m = LineMask(lam0=lam0, note=note)
        m.set_fit(lmin, lmax, mset.pad)
        mset.masks.append(m)
        self.active = len(mset.masks) - 1
        self._changed()
        return self.active

    def remove(self, index):
        if 0 <= index < len(self.masks):
            self.masks.pop(index)
            if not self.masks:
                self.active = None
            elif self.active is not None and self.active >= len(self.masks):
                self.active = len(self.masks) - 1
            self._changed()

    def set_active(self, index):
        self.active = index if (index is not None and 0 <= index < len(self.masks)) else None
        self._changed()

    def set_fit(self, index, lmin, lmax):
        """Resize mask ``index``'s fit window; its synth window follows the set ``pad``."""
        if 0 <= index < len(self.masks):
            self.masks[index].set_fit(lmin, lmax, self.mset.pad)
            self._changed()

    def set_element(self, element):
        """Set which element the current set's windows diagnose (persisted to its file)."""
        value = str(element or "").strip()
        self.mset.element = value.capitalize() if value else ""
        self._changed()

    def set_pad(self, pad):
        """Re-derive every mask's synth/visual window from a new padding (Å)."""
        self.mset.repad(pad)
        self._changed()

    def set_renorm(self, index, on):
        if 0 <= index < len(self.masks):
            self.masks[index].renorm = bool(on)
            self._changed()

    # --- session persistence ------------------------------------------------------------
    def snapshot(self) -> dict:
        """The ``masks`` section of a workspace slot: every set (baseline included, so the
        clean/modified state survives), the current set and the active mask index."""
        return {"sets": [s.to_dict() for s in self.sets],
                "current": self.current_uid, "active": self.active}

    def restore(self, sets, current=None, active=None):
        """Replace every set with the prepared ``sets`` (``LineMaskSet`` objects carrying
        their saved uids; see ``session.prepare``). An empty list leaves one fresh scratch
        set, as a new document has. ``current`` must be a uid of ``sets`` or ``None`` (the
        last set); ``active`` an index into the current set's masks or ``None``."""
        self.sets = list(sets) or [LineMaskSet()]
        for s in self.sets:
            if not s.uid:
                s.uid = next(self._uids)
        self._uids = itertools.count(max(s.uid for s in self.sets) + 1)
        self.current_uid = current if self.find(current) is not None else self.sets[-1].uid
        self.active = (active if active is not None and 0 <= active < len(self.masks)
                       else None)
        self._changed()
