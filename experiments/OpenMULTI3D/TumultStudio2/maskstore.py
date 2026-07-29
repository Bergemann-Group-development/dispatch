"""The shared LineMask set for TumultStudio2.

One :class:`MaskStore` holds the app's current :class:`linemask.LineMaskSet` — the single
source of truth every consumer reads: the canvas draws it, the LineMasks tool edits it, and
(later) the abundance / synthesize-fit workflows adopt it. Mutations broadcast the context's
``"masks-changed"`` event so every subscriber re-syncs, exactly like the existing
``"changed"`` / ``"filter-state"`` events.

The store is deliberately GUI-free (it only touches the context's event bus) so it can be
unit-tested and serialised into a session snapshot on its own.
"""
from __future__ import annotations

from linemask import DEFAULT_PAD, LineMask, LineMaskSet, parse_linemask_file


class MaskStore:
    def __init__(self, ctx):
        self.ctx = ctx
        self.mset = LineMaskSet()
        self.active = None      # index of the mask being edited on the canvas, or None

    # --- queries ------------------------------------------------------------------------
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
        """Replace the set with the masks parsed from ``path`` (windows only, no measuring)."""
        self.mset = parse_linemask_file(path, pad=self.mset.pad)
        self.active = 0 if self.masks else None
        self._changed()
        return self.mset

    def clear(self):
        self.mset = LineMaskSet(pad=self.mset.pad)
        self.active = None
        self._changed()

    def add(self, lmin, lmax, *, lam0=float("nan"), note="") -> int:
        """Append a new mask over ``[lmin, lmax]`` (synth window derived from ``pad``); select it."""
        m = LineMask(lam0=lam0, note=note)
        m.set_fit(lmin, lmax, self.mset.pad)
        self.masks.append(m)
        self.active = len(self.masks) - 1
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

    def set_pad(self, pad):
        """Re-derive every mask's synth/visual window from a new padding (Å)."""
        self.mset.repad(pad)
        self._changed()

    def set_renorm(self, index, on):
        if 0 <= index < len(self.masks):
            self.masks[index].renorm = bool(on)
            self._changed()

    def set_cont_pts(self, index, pts):
        if 0 <= index < len(self.masks):
            self.masks[index].cont_pts = tuple(tuple(float(c) for c in p) for p in pts)
            self._changed()

    # --- session persistence ------------------------------------------------------------
    def snapshot(self) -> dict:
        return {"mset": self.mset.to_dict(), "active": self.active}

    def restore(self, data):
        data = data or {}
        self.mset = LineMaskSet.from_dict(data.get("mset")) if data.get("mset") else LineMaskSet()
        act = data.get("active")
        self.active = act if (act is not None and 0 <= act < len(self.masks)) else None
        self._changed()
