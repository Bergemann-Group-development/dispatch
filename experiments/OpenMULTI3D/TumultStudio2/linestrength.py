"""The strength ladder: which source answers for an object and a window.

The strengths themselves, and the reference atmospheres they are computed on, are
:mod:`tumult.linestrength` and :mod:`tumult.refstars`. What is here is the choice between
sources and the polling shell the UI needs to make that choice without blocking a document:

1. **Computed** — HOP on a reference atmosphere, with the configured line list. Which
   atmosphere is the line-list panel's own dropdown: a named class, or ``auto``, which
   follows the spectrum's own (Teff, log g, [Fe/H]) — from a synthesis's parameters or a run
   folder's namelist — to the nearest class, and falls back to the G dwarf for a spectrum
   that says nothing.
2. **The master line list and its proxy** — when HOP cannot be reached (no ``libhop``) or
   the band has not been computed yet. The only rung whose ordering is not physical, and it
   is labelled ``approximate`` wherever it is shown.

Both hand back the same frame — ``ion``/``lam0``/``evi``/``loggf``/``strength`` ascending in
lam0 — so callers never branch on which rung answered.
"""
from __future__ import annotations

import logging
from pathlib import Path

from tumult import refstars
from tumult.linelist import MASTER_PATH, load_master
from tumult.linestrength import (DEFAULT_FEH, DEFAULT_LOGG, band_bounds, band_ready,
                                 cached_band, compute_band, run_params)
# the shared "columns every rung returns" normalisation: the proxy rung has to hand back
# the same frame shape as a computed band, so both go through the one implementation
from tumult.linestrength import _sorted_frame

class Source:
    """One rung's answer: the lines, what to call it, and whether it can be trusted to rank.

    ``key`` is what the frame *is*, as opposed to ``label``, which is what a footer calls it.
    A consumer caching the frame has to key on this: the label names the star but not the
    band, so two answers a pan apart look identical by it while holding different lines.
    """

    __slots__ = ("lines", "label", "kind", "params", "key")

    def __init__(self, lines, label, kind, params=None, key=None):
        self.lines = lines
        self.label = label
        self.kind = kind       # "computed" | "proxy"
        self.params = params
        self.key = key if key is not None else (kind, label)

    @property
    def approximate(self):
        return self.kind == "proxy"

    def __repr__(self):
        n = 0 if self.lines is None else len(self.lines)
        return f"<Source {self.kind} {self.label!r} n={n}>"


class LineStrengths:
    """The strength ladder, shared by the annotations and the line-list panel."""

    def __init__(self, master_path=MASTER_PATH, linelist=None,
                 refstar=refstars.AUTO_KEY):
        self.master_path = str(master_path)
        # what HOP computes strengths *for*: the master list unless told otherwise, so the
        # ranking and the table describe the same transitions
        self.linelist = str(linelist or master_path)
        # the chosen reference atmosphere, or AUTO_KEY to follow each spectrum's own. Via
        # the setter, because this normally arrives from a stored preference and a key that
        # has since left the table has to land on auto, not on an unrelated star.
        self.refstar = refstars.AUTO_KEY
        self.set_refstar(refstar)
        self._proxy = None

    # --- configuration ---------------------------------------------------------------
    def set_master_path(self, path):
        self.master_path = str(path)
        self.linelist = str(path)
        self._proxy = None

    def set_refstar(self, key):
        """Choose the reference atmosphere, by :mod:`refstars` key or ``auto``.

        Returns what stuck, so a caller can persist exactly what it set. An unknown key
        falls back to ``auto`` rather than to some arbitrary star: a stale preference should
        follow the spectrum, not silently rank it as something it is not.
        """
        key = str(key)
        self.refstar = key if (key == refstars.AUTO_KEY or key in refstars.BY_KEY) \
            else refstars.AUTO_KEY
        return self.refstar

    # --- what star is this? ------------------------------------------------------------
    def params_of(self, obj):
        """``(teff, logg, feh)`` stated by ``obj`` — its synthesis parameters or its run's
        namelist — or ``None`` when it says nothing about its star."""
        meta = getattr(obj, "meta", None) or {}
        stellar = meta.get("stellar") if isinstance(meta, dict) else None
        if isinstance(stellar, dict) and stellar.get("teff"):
            try:
                return (float(stellar["teff"]),
                        float(stellar.get("logg", DEFAULT_LOGG)),
                        float(stellar.get("feh", DEFAULT_FEH)))
            except (TypeError, ValueError):
                pass
        source = getattr(obj, "source", None)
        if source and len(source) > 1 and source[0] == "run" and source[1]:
            found = run_params(source[1])
            if found:
                return found
        return None

    def star_for(self, obj):
        """The :class:`refstars.RefStar` the strengths for ``obj`` are computed on.

        A chosen class wins outright — that is what choosing one is for, and it lets a
        spectrum be ranked deliberately against something other than itself. Under ``auto``
        the spectrum's own parameters pick the nearest class, and a spectrum that states none
        gets the solar one, which is the honest default for an unlabelled observation.
        """
        if self.refstar != refstars.AUTO_KEY:
            return refstars.get(self.refstar)
        params = self.params_of(obj)
        return refstars.nearest(*params) if params else refstars.get(refstars.DEFAULT_KEY)

    # --- the ladder ------------------------------------------------------------------
    def lines_for(self, obj, lo, hi, allow_compute=False):
        """The best :class:`Source` for ``obj`` over ``[lo, hi]`` Å.

        With ``allow_compute`` false — what a repaint on the document thread passes — a band
        that has not been computed is not waited for; the proxy answers instead and
        :meth:`warm` does the computing in the background, after which the next repaint gets
        the real thing.
        """
        star = self.star_for(obj)
        band = band_bounds(lo, hi)
        if allow_compute:
            frame = compute_band(star, self.linelist, *band)
        else:
            # Never enter compute_band from a repaint, even if another thread
            # evicts a previously ready band between the readiness check and read.
            frame = cached_band(star, self.linelist, *band)
        if frame is not None and len(frame):
            return Source(frame, f"HOP · {star.label} ({star.summary})",
                          "computed", (star.teff, star.logg, star.feh),
                          key=("hop", star.key, self.linelist, band))
        return self._proxy_source()

    def warm(self, obj, lo, hi):
        """Compute what :meth:`lines_for` wants, blocking. Call off the document thread."""
        return self.lines_for(obj, lo, hi, allow_compute=True)

    def pending(self, obj, lo, hi):
        """True when :meth:`warm` would change the answer — i.e. the band is not computed."""
        return not band_ready(self.star_for(obj), self.linelist, *band_bounds(lo, hi))

    def describe(self, obj, lo, hi):
        """One short line naming the source in play, for the panel footer."""
        src = self.lines_for(obj, lo, hi)
        return ("strengths: approximate (computing…)" if src.approximate
                else f"strengths: {src.label}")

    # --- the proxy rung ------------------------------------------------------------------
    def _proxy_source(self):
        """The configured master list, with whatever strength column it has.

        Normally the plain 6-column form, so this is the Boltzmann proxy ``load_master``
        falls back to — the one ordering that does not survive a change of species. A list
        that happens to carry real strengths (a ``linelist_out.txt``) keeps them, and then
        the source stops being approximate.
        """
        if self._proxy is None:
            try:
                master = load_master(self.master_path)
            except Exception:
                logging.exception("Line strengths: could not read %s", self.master_path)
                self._proxy = Source(None, "unavailable", "proxy")
            else:
                proxy = bool(master.attrs.get("strength_is_proxy", True))
                name = Path(self.master_path).name
                self._proxy = Source(_sorted_frame(master),
                                     f"{name} (approximate)" if proxy else name,
                                     "proxy" if proxy else "computed",
                                     key=("master", self.master_path))
        return self._proxy
