"""The reference stars: named sets of stellar parameters, shared by the GUI.

What a reference star is
------------------------
One entry is a label and ``(Teff, log g, [Fe/H], vmic)``. The list serves three places:
the *Reference star* dropdown of every synthesis setup, which fills the stellar-parameter
fields from an entry; the ionisation-balance entry point, which starts a solve from one;
and the line-strength ranking, which computes strengths on the reference atmosphere and
picks the entry nearest to a spectrum's own parameters when none is chosen.

The list is editable (Preferences → Reference stars): the GUI stores it and installs it
here with :func:`set_stars`, so every consumer reads one catalogue. :data:`DEFAULT_STARS`
is what a fresh configuration starts from — nine classes spanning 2500-6500 K and log g
0-5, the cool half of the HR diagram where a line list has anything to say.

Where the atmosphere comes from
-------------------------------
An entry carries parameters only. The T/rho structure a strength calculation needs is looked
up on the registered MARCS grid by :func:`tumult.linestrength.atmosphere_for`: the node
nearest to the entry in the same scaled distance :func:`nearest` uses, read at the eight
depths :func:`tumult.linestrength.read_marcs` samples. Adding a star is therefore just
naming its parameters.

Choosing between them
---------------------
:func:`nearest` picks by scaled distance: 1000 K of Teff, 1 dex of gravity, 1 dex of
metallicity and 1 km/s of microturbulence weighted alike, so a spectrum that states its own
parameters lands on a sensible entry without anyone choosing one. The names are there for
the spectra that state nothing, which is most observed ones: a class is something you can
pick off a dropdown from what you know about the star.
"""
from __future__ import annotations

import re

import numpy as np

PARAM_KEYS = ("teff", "logg", "feh", "vmic")


class RefStar:
    """One reference star: what to call it and its ``(teff, logg, feh, vmic)``."""

    __slots__ = ("key", "label", "teff", "logg", "feh", "vmic")

    def __init__(self, key, label, teff, logg, feh, vmic=1.0):
        self.key = str(key)
        self.label = str(label)
        self.teff = float(teff)
        self.logg = float(logg)
        self.feh = float(feh)
        self.vmic = float(vmic)

    @property
    def params(self):
        """``{teff, logg, feh, vmic}`` as floats, the order of :data:`PARAM_KEYS`."""
        return {k: getattr(self, k) for k in PARAM_KEYS}

    @property
    def summary(self):
        """``Teff 5750, log g 4.50, [Fe/H] +0.00`` — the parameters, for a footer."""
        return f"Teff {self.teff:.0f}, log g {self.logg:.2f}, [Fe/H] {self.feh:+.2f}"

    def matches(self, teff, logg, feh, vmic):
        """Whether these values are this star's, to the precision the fields show."""
        return (abs(self.teff - float(teff)) < 0.5 and abs(self.logg - float(logg)) < 0.005
                and abs(self.feh - float(feh)) < 0.005 and abs(self.vmic - float(vmic)) < 0.05)

    def as_dict(self):
        return {"key": self.key, "label": self.label, **self.params}

    @classmethod
    def from_dict(cls, d):
        return cls(d["key"], d.get("label") or d["key"], d["teff"], d["logg"], d["feh"],
                   d.get("vmic", 1.0))

    def __repr__(self):
        return f"<RefStar {self.key} {self.summary}, vmic {self.vmic:.1f}>"


# the classes a fresh configuration starts from; vmic follows the MARCS convention of
# 1 km/s for dwarfs, 2 for giants and 5 for supergiants
DEFAULT_STARS = [
    RefStar("gdwarf", "G dwarf (solar)", 5750, 4.5, +0.00, 1.0),
    RefStar("fdwarf", "F dwarf", 6500, 4.0, +0.00, 1.0),
    RefStar("kdwarf", "K dwarf", 4500, 4.5, +0.00, 1.0),
    RefStar("mdwarf", "M dwarf", 3500, 5.0, +0.00, 1.0),
    RefStar("turnoff", "Turn-off (metal-poor)", 6250, 4.0, -2.00, 1.0),
    RefStar("clump", "Red clump", 4750, 2.5, +0.00, 2.0),
    RefStar("rgb", "Red giant branch", 4500, 2.0, +0.00, 2.0),
    RefStar("rgbmp", "Metal-poor giant", 4500, 1.5, -2.00, 2.0),
    RefStar("rsg", "Red supergiant", 3800, 0.0, +0.00, 5.0),
]

DEFAULT_KEY = "gdwarf"
AUTO_KEY = "auto"                # follow the spectrum's own parameters, when it has any
CUSTOM_KEY = "custom"            # a setup's own numbers, matching no entry

_stars = list(DEFAULT_STARS)


def stars():
    """The current catalogue, in display order."""
    return list(_stars)


def by_key():
    return {s.key: s for s in _stars}


def set_stars(entries):
    """Install the catalogue: :class:`RefStar` objects or their ``as_dict`` mappings.
    Entries repeating an earlier key are dropped."""
    global _stars
    seen, out = set(), []
    for e in entries:
        star = e if isinstance(e, RefStar) else RefStar.from_dict(e)
        if star.key in seen:
            continue
        seen.add(star.key)
        out.append(star)
    _stars = out


def slug(label, existing=()):
    """A key for a new entry named ``label``: lower-case ASCII words, made unique against
    ``existing`` with a numeric suffix."""
    base = re.sub(r"[^a-z0-9]+", "_", str(label).lower()).strip("_") or "star"
    key, n = base, 1
    while key in existing or key in (AUTO_KEY, CUSTOM_KEY):
        n += 1
        key = f"{base}{n}"
    return key


def options(auto_label="Auto (from spectrum)"):
    """``[(key, label)]`` for a select field, ``auto`` first."""
    return [(AUTO_KEY, auto_label)] + [(s.key, s.label) for s in _stars]


def get(key):
    """The star named ``key``, else the default entry (or the first one when the default
    was removed), else ``None`` on an empty catalogue. ``auto`` is not a star — it
    resolves against a spectrum's parameters, so callers wanting one go through
    :meth:`linestrength.LineStrengths.star_for`."""
    table = by_key()
    if str(key) in table:
        return table[str(key)]
    if DEFAULT_KEY in table:
        return table[DEFAULT_KEY]
    return _stars[0] if _stars else None


def matching(teff, logg, feh, vmic):
    """The entry whose parameters these are (see :meth:`RefStar.matches`), else ``None``."""
    return next((s for s in _stars if s.matches(teff, logg, feh, vmic)), None)


def nearest(teff, logg, feh, vmic=None):
    """The entry closest to ``(teff, logg, feh[, vmic])``, or ``None`` on an empty list.

    Distances are scaled so the axes are comparable: 1000 K of Teff, 1 dex of gravity,
    1 dex of metallicity and 1 km/s of microturbulence count the same. Teff dominates in
    practice, which is what lets a nine-entry table stand in for a grid. ``vmic`` is
    ignored when not given, as by a spectrum stating only three parameters.
    """
    if not _stars:
        return None
    t = np.array([s.teff for s in _stars])
    g = np.array([s.logg for s in _stars])
    z = np.array([s.feh for s in _stars])
    d = ((t - float(teff)) / 1000.0) ** 2 + (g - float(logg)) ** 2 + (z - float(feh)) ** 2
    if vmic is not None:
        d = d + (np.array([s.vmic for s in _stars]) - float(vmic)) ** 2
    return _stars[int(np.argmin(d))]
