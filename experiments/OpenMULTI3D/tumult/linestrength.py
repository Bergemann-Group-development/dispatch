"""Where a line's *strength* comes from — the one number that decides which transitions
are worth naming and which head a line-list table.

The problem this solves
-----------------------
The obvious proxy, ``lam0 · 10^loggf · exp(-E_low·EK/T)`` at some fixed T, is not a strength
at all: it carries no element abundance, no ionization fraction, no molecular dissociation
and no partition function, so it cannot be compared between two species. Ranked by it, a
two-Angstrom window of the solar spectrum is headed by TiO and C2 bands that are not in the
Sun at all, and its ordering correlates with the true one at Spearman rho ~ 0.2 — no better
than chance.

HOP computes the right quantity per line while it builds the opacity::

    ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))

— line-centre extinction over continuum extinction, maximised over the atmosphere. It is
dimensionless, it compares across species, and it is the criterion the RT code itself uses
to decide which lines are worth computing.

Computed, not stored
--------------------
That number used to be readable only from the ``linelist_out.txt`` a full run writes, which
meant keeping a file per reference star — up to a gigabyte for a cool one, and silently
wrong the moment the configured line list changed, since a dump only ever describes the list
its run used. ``hop_linelist_count`` / ``hop_get_linelist`` / ``hop_linelist_species`` in the
standalone C API (``HOP/standalone/hop_c_api.f90``) hand the same numbers back from a solved
state in memory. :func:`compute_band` is one EOS solve on a reference atmosphere plus one
opacity call over a band of spectrum, and it returns the strengths of every line in it.

The atmosphere is one of a handful of named classes held as a table of eight (T, rho) pairs
in :mod:`tumult.refstars`, not read off the MARCS grid. Eight depths is all a strength needs
(see :func:`read_marcs`), which makes the whole atmosphere sixteen floats: too little to
justify finding, name-parsing and opening one of fifteen thousand files.

Measured against ``data/sed_3000-9000``'s own dump, over 5000-5002 A with that run's line
list: all 155 lines matched, **Spearman 0.978**, and the ratio to the run's values is a near
constant 1.21 (10-90% spread 1.21-1.33) — an offset from using a 1D MARCS structure where
the run had its own, which a ranking does not care about.

Cost, for the single-file master list: 0.06 s for a 20 A band, 0.08 s for 100 A, 0.32 s for
500 A, and 3.3 s for the whole 4000-9000 A at once, because only :data:`NDEPTH` depths are
passed. A folder of per-species lists is ~5x that, dominated by reading the files. Bands are
snapped to a grid and cached process-wide, so panning inside one costs nothing at all.

Every band returns the same frame — ``ion``/``lam0``/``evi``/``loggf``/``strength``
ascending in lam0 — as :func:`tumult.linelist.load_master` does with its proxy column, so a
caller that has to fall back never branches on which of the two answered.

Threading
---------
HOP keeps process-wide global state and refuses a second live instance. The band cache is
module-level and the native call is serialised on
:data:`tumult.synthesis.runtime.HOP_LOCK`, shared by every caller in the process: two
sessions looking at the same star reuse one computation rather than fighting over the
library.
"""
from __future__ import annotations

import json
import logging
import re
import threading
from pathlib import Path

import numpy as np

from tumult.synthesis.api import PROVENANCE_NAME
from tumult.synthesis.runtime import HOP_LOCK

DEFAULT_LOGG = 4.44          # filled in when a spectrum states a Teff but not these
DEFAULT_FEH = 0.0
BAND_AA = 100.0              # bands snap to this grid, so a pan usually stays inside one
BAND_STEPS = 16              # a view spans at most this many grid steps; the step doubles
                             # until it does, so wide views share a band instead of each
                             # asking for its own (see band_bounds)
MAX_BAND_AA = 6000.0         # widest band computed in one go (cost is linear in width) —
                             # wide enough that a whole optical spectrum is one band
GRID_AA = 0.005              # sampling: must resolve a line core, or peaks are under-read
DLAM_AA = 2.0                # wing search radius handed to HOP's line-list reader
MAX_BANDS = 8                # bands kept in the process-wide cache
NDEPTH = 8                   # depths of the reference atmospheres (see read_marcs)
PARAMS_NAME = "params.nml"

_ATMOS_RE = re.compile(r"ATMOS_FILE='\s*([^']*?)\s*'", re.I)
_FEH_RE = re.compile(r"\bFEH=\s*([-+\d.Ee]+)")

# The bands are identical for every caller, so the cache is process-wide; the native call
# itself is serialised on the shared HOP_LOCK. A separate, short-lived cache lock lets a
# caller check readiness while another thread is inside HOP.
_bands_lock = threading.Lock()
_bands = {}      # (star, linelist, lo, hi) -> frame, most recent last


def _sorted_frame(df):
    """The columns every rung returns, ascending in lam0, with a usable index."""
    cols = [c for c in ("ion", "lam0", "evi", "loggf", "strength") if c in df.columns]
    out = df[cols].copy()
    out["ion"] = out.ion.astype(str).str.strip()
    return out.sort_values("lam0").reset_index(drop=True)


# --- the model atmosphere -----------------------------------------------------------------
def read_marcs(path):
    """``(T[K], rho[g/cm3])`` down a MARCS model — the two columns the EOS needs.

    Nothing in the app calls this: the atmospheres it uses are the table in :mod:`refstars`.
    This is where that table came from and how it is extended — point it at a MARCS model and
    paste the two arrays into a new :class:`refstars.RefStar`.

    The structure is in two tables keyed by the same depth index: T in the first, density in
    the second. A strength is a maximum over the cells of the state, so what the sampling has
    to preserve is *where that maximum lives* — and it lives anywhere: high-excitation lines
    peak in the deep hot layers, strong resonance and molecular ones at the very top. Only
    :data:`NDEPTH` layers are passed, spread over the model's whole index range, which is
    what makes this cheap; the cost of a band is roughly linear in the count, and the full
    56 is about 4x this for nothing. Measured over 5000-5100 Angstrom on a solar model
    against all 56 layers: Spearman 1.0000, every line matched, strengths within 10%.

    Spanning the range is the part that matters, not the count. Four layers over
    lgTau5 = -3..0 — an obvious-looking choice, and a plausible one for a *representative*
    depth — is wrong here for exactly this reason: it clips both ends of the maximum, so the
    median line comes out 2.7x too weak, the worst 30x (1000x on a 3500 K model), and a tenth
    of the lines fall under HOP's own cutoff and are not returned at all.
    """
    lines = Path(path).read_text(errors="ignore").split("\n")
    n = int(next(l for l in lines if "Number of depth points" in l).split()[0])
    i1 = next(i for i, l in enumerate(lines) if "lgTauR" in l and "Pe" in l)
    i2 = next(i for i, l in enumerate(lines) if "lgTauR" in l and "Density" in l)
    keep = np.unique(np.linspace(0, n - 1, min(NDEPTH, n)).round().astype(int))
    T = np.array([float(lines[i1 + 1 + k].split()[4]) for k in keep])
    rho = np.array([float(lines[i2 + 1 + k].split()[3]) for k in keep])
    return T, rho


# --- the computation ----------------------------------------------------------------------
def band_bounds(lo, hi):
    """``(lo, hi)`` of the cache band covering ``[lo, hi]``, snapped to :data:`BAND_AA`.

    Snapping is what makes panning free: neighbouring views ask for the same band and hit the
    cache. The step is :data:`BAND_AA` while the view spans at most :data:`BAND_STEPS` of
    them, and doubles until it does — so a narrow view keeps the fine grid it always had,
    while every whole-spectrum view lands on the same coarse band and computes once.

    A centred crop to :data:`MAX_BAND_AA` used to apply from 1000 Å up, on the reasoning that
    a frame holds only a few tens of names anyway. That is true, but they belong *spread
    across the view*: cropping to the middle left the annotation overlay drawing every name
    it had in a narrow clump at the centre of a full-spectrum plot, with the rest of the
    range bare. The cap is now wide enough to cover an optical spectrum end to end and only
    bites on something far wider.
    """
    lo, hi = float(min(lo, hi)), float(max(lo, hi))
    if hi - lo > MAX_BAND_AA:
        mid = 0.5 * (lo + hi)
        lo, hi = mid - MAX_BAND_AA / 2, mid + MAX_BAND_AA / 2
    step = BAND_AA
    while hi - lo > BAND_STEPS * step:
        step *= 2
    return np.floor(lo / step) * step, np.ceil(hi / step) * step


def compute_band(star, linelist, lo, hi):
    """Strengths of every line in ``[lo, hi]`` Å, from HOP on the ``star`` atmosphere.

    ``star`` is a :class:`tumult.refstars.RefStar` — it carries its own T and rho, so this touches
    no files but the line list. Blocking, and seconds at worst for a wide band, so call it
    off the document thread. Serialised on :data:`tumult.synthesis.runtime.HOP_LOCK` because HOP is a
    process-wide singleton. A separate short-lived cache lock lets the UI check readiness while
    another thread computes a band.
    """
    key = (star.key, str(linelist), lo, hi)
    with _bands_lock:
        if key in _bands:
            return _bands[key]
    with HOP_LOCK:
        with _bands_lock:
            if key in _bands:
                return _bands[key]
        try:
            import tumult as tm
            from tumult.physics import air2vac
            lam = air2vac(np.arange(lo, hi, GRID_AA))
            hop = tm.HopEOS(feh=0.0)
            try:
                hop.load_state(star.T, rho=star.rho)
                df = hop.line_strengths(lam, str(linelist), dlam=DLAM_AA, vmic=1.0)
            finally:
                hop.close()
            frame = _sorted_frame(df[(df.lam0 >= lo) & (df.lam0 <= hi)])
        except Exception:
            logging.exception("Line strengths: HOP failed for %s over %.1f-%.1f",
                              star.key, lo, hi)
            frame = None
        with _bands_lock:
            _bands[key] = frame
            while len(_bands) > MAX_BANDS:
                _bands.pop(next(iter(_bands)))
        return frame


def band_ready(star, linelist, lo, hi):
    """Whether :func:`compute_band` would return without computing anything."""
    with _bands_lock:
        return (star.key, str(linelist), lo, hi) in _bands


def run_params(folder):
    """``(teff, logg, feh)`` of a run folder, for picking the reference atmosphere.

    A run the synthesis API launched states them in its own ``provenance.json``. Any other
    run folder -- one produced by hand, outside the cache -- has only its ``params.nml``,
    whose atmosphere file name carries the real values; the folder name is rounded (``t58``
    for 5750 K) and is not trusted.
    """
    try:
        record = json.loads((Path(folder) / PROVENANCE_NAME).read_text())
        star = record["star"]
        return (float(star["teff"]), float(star["logg"]), float(star["feh"]))
    except (OSError, ValueError, KeyError, TypeError):
        pass
    try:
        txt = (Path(folder) / PARAMS_NAME).read_text(errors="ignore")
    except OSError:
        return None
    m = _ATMOS_RE.search(txt)
    if not m:
        return None
    a = re.search(r"[ps](\d+)_g([+-]?[\d.]+)", Path(m.group(1)).name)
    if not a:
        return None
    feh = _FEH_RE.search(txt)
    return (float(a.group(1)), float(a.group(2)),
            float(feh.group(1)) if feh else DEFAULT_FEH)
