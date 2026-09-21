"""Names of the objects a calculation adds to the document.

A generated name is fixed when the calculation is submitted: the object it delivers, the
workflow record's label and the Processes row all carry it, so a row says which object it
will update. The name states what tells alternatives apart (the target, the window, the
atmosphere, LTE or NLTE, abundance overrides); the Info tab and the result panels hold the
full provenance. :func:`unique` keeps names distinct across the document with a ``#n``
suffix; names already in a loaded workspace are kept as saved.

Examples::

    T5750 g4.40 m+0.00 ξ1.00 · 5000–5001 Å · MARCS · LTE
    obs.fits fit · fe-lines · MARCS · NLTE fe.atom · [Mg/Fe]+0.20
    Fe abundance · obs.fits · fe-lines · MARCS · LTE
    Fe corrections · T5750 g4.40 m+0.00 ξ1.00 · MARCS
"""
from __future__ import annotations

import re
from pathlib import Path

from tumult.departures import DepartureGrid, DepartureRequest, Departures

SEPARATOR = " · "
_SUFFIX = re.compile(r"^(.*) #(\d+)$")


def unique(name, taken):
    """``name``, or ``name #n`` with the smallest ``n >= 2`` not in ``taken``."""
    taken = set(taken)
    if name not in taken:
        return name
    match = _SUFFIX.match(name)
    stem = match.group(1) if match else name
    n = 2
    while f"{stem} #{n}" in taken:
        n += 1
    return f"{stem} #{n}"


def _number(value):
    """``value`` with at most two decimals and no trailing zeros."""
    return f"{float(value):.2f}".rstrip("0").rstrip(".")


def star(teff, logg, feh, vmic=None):
    text = f"T{float(teff):.0f} g{float(logg):.2f} m{float(feh):+.2f}"
    return text if vmic is None else f"{text} ξ{float(vmic):.2f}"


def span(windows):
    """``lo–hi Å`` over every ``(lo, hi)`` window."""
    lo = min(float(w[0]) for w in windows)
    hi = max(float(w[1]) for w in windows)
    return f"{_number(lo)}–{_number(hi)} Å"


def atmosphere(request, grid=None):
    """The picked atmosphere file, else the interpolated grid's registered name."""
    path = getattr(request.atmosphere, "path", None)
    if path:
        return Path(str(path)).name
    grid = grid or {}
    return str(grid.get("name") or grid.get("id") or request.atmosphere.format)


def atoms(request):
    """The model atom files of the request's departures, as strings (none for LTE)."""
    departures = request.departures
    if isinstance(departures, DepartureRequest):
        return [str(atom) for atom in departures.atoms]
    if isinstance(departures, (Departures, DepartureGrid)):
        return list(departures.atom_settings[0])
    return []


def mode(request):
    """``LTE``, or ``NLTE`` with the model atoms' file names."""
    if not request.is_nlte:
        return "LTE"
    names = ",".join(Path(atom).name for atom in atoms(request))
    return f"NLTE {names}" if names else "NLTE"


def abundances(request, skip=()):
    """The request's abundance overrides, except the elements in ``skip``."""
    values = request.abundances
    if not values:
        return ""
    relative = values.relative
    parts = [(f"[{el}/Fe]{value:+.2f}" if relative else f"A({el}){_number(value)}")
             for el, value in values.items() if el not in skip]
    return " ".join(parts)


def join(*parts):
    return SEPARATOR.join(str(p) for p in parts if p)


def synthesis(request, grid=None, window=None):
    """A synthesis: the star first, then the window (``window`` names a line mask)."""
    return join(star(*request.star.tuple()), window or span(request.windows),
                atmosphere(request, grid), mode(request), abundances(request))


def fit(target, request, grid=None, window=None, fitted=()):
    """A spectrum fit of ``target``; abundances being fitted are not overrides."""
    return join(f"{target} fit", window or span(request.windows),
                atmosphere(request, grid), mode(request), abundances(request, fitted))


def analysis(element, kind, target, mask, request, grid=None):
    """A line-mask analysis (``kind`` is ``abundance`` or ``balance``) of ``target``."""
    return join(f"{element} {kind}", target, mask, atmosphere(request, grid), mode(request))


def corrections(element, request, grid=None):
    """NLTE corrections of ``element`` at the request's star."""
    return join(f"{element} corrections", star(*request.star.tuple()),
                atmosphere(request, grid))


def departures(element, request, grid=None):
    """A departure solve of ``element`` on the request's atmosphere."""
    return join(f"{element} departures", star(*request.star.tuple()),
                atmosphere(request, grid))
