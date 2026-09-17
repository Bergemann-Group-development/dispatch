"""Placing a continuum under an observed spectrum, and dividing by it.

A synthetic spectrum comes out of the code already normalised; an observation does not, and
the abundance and equivalent-width machinery downstream assumes a flux of roughly unity
between the lines. Everything here answers one question — where is the continuum — in the
three ways the analysis needs it:

* :func:`flat_level` — one robust percentile of the flux. The default when nothing better is
  known. A high percentile rather than the maximum, so a cosmic ray or a hot pixel does not
  set the level for the whole spectrum.
* :func:`merge_knots` / :func:`evaluate` — a monotonic spline through hand-placed knots,
  optionally overridden inside local windows by their own colinear knots. Monotone-cubic
  (pchip) rather than a plain cubic because a continuum must not overshoot between knots:
  an overshoot puts the flux above the continuum and turns a line into emission.
* :func:`reject_outliers` — the preprocessing a spline is fitted *to*: low-pass the spectrum
  to the resolution of interest, then drop points sitting far above a running mean, which is
  where absorption lines and noise spikes live.

These are pure array operations over ``(x, y)``; which spectrum they apply to, and how the
knots got there, belongs to the caller.
"""
from __future__ import annotations

import numpy as np

from .spectrum import Spectrum
from .utils import pchip_mono

# A high percentile of the flux, not its maximum: robust to cosmic rays and hot pixels while
# still sitting in the line-free part of a normal stellar spectrum.
NORM_PERCENTILE = 98.0

# Local windows are narrower, so proportionally more of a window is line rather than
# continuum: the cut that isolates the continuum there has to sit lower than the global one.
SEED_PERCENTILE = 90.0


def flat_level(y, percentile=NORM_PERCENTILE):
    """The flat continuum level for flux ``y``: a robust percentile, NaNs ignored."""
    return float(np.nanpercentile(np.asarray(y, float), percentile))


def flat_continuum(y, percentile=NORM_PERCENTILE):
    """:func:`flat_level` broadcast to ``y``'s shape."""
    y = np.asarray(y, float)
    return np.full(y.shape, flat_level(y, percentile))


def upper_envelope_fit(x, y, percentile=SEED_PERCENTILE):
    """A straight line through the upper-``percentile`` points of ``(x, y)``.

    The continuum sits above the absorption, so the points worth fitting are the high ones.
    Returns ``(slope, intercept)``, or ``None`` when there is nothing to fit — fewer than two
    points survive the cut, or they all sit at one x, which no line can be fitted through.
    A caller that gets ``None`` should fall back to :func:`flat_level`.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 4:
        return None
    top = y >= np.nanpercentile(y, percentile)
    if top.sum() < 2 or np.ptp(x[top]) <= 0:
        return None
    slope, intercept = np.polyfit(x[top], y[top], 1)
    return float(slope), float(intercept)


def merge_knots(x, y, base_knots=None, renorm_windows=(), percentile=NORM_PERCENTILE):
    """The spline knots for a continuum: ``base_knots`` merged with local-window knots.

    ``base_knots`` is the global spline as ``(kx, ky)``, or ``None`` for a flat percentile
    level spanning ``x``. Each entry of ``renorm_windows`` is ``(smin, smax, knots)``, whose
    ``knots`` are ``(x, y)`` pairs governing that window. Base knots strictly inside a window
    are dropped so the local continuum governs there, while the padding beyond blends into
    the rest of the spline.

    Returns ``(kx, ky)`` sorted and strictly increasing in ``kx`` — :func:`pchip_mono`
    requires that, and colinear window knots readily collide.
    """
    x = np.asarray(x, float)
    if base_knots is not None:
        pts = list(zip((float(v) for v in base_knots[0]),
                       (float(v) for v in base_knots[1])))
    else:
        lvl = flat_level(y, percentile)
        pts = [(float(np.min(x)), lvl), (float(np.max(x)), lvl)]
    for smin, smax, knots in (renorm_windows or []):
        if not knots:
            continue
        pts = [(px, py) for px, py in pts if not (smin < px < smax)]
        pts += [(float(a), float(b)) for a, b in knots]
    pts.sort(key=lambda p: p[0])
    kx, ky = [], []
    for px, py in pts:
        if kx and px <= kx[-1]:
            px = kx[-1] + 1e-6      # pchip_mono needs strictly-increasing x
        kx.append(px)
        ky.append(py)
    return np.asarray(kx), np.asarray(ky)


def evaluate(kx, ky, x):
    """The monotonic spline through ``(kx, ky)``, sampled at ``x``.

    ``x`` is clamped to the knot span so the continuum holds flat past the outermost knots
    rather than diverging on the cubic extrapolation — an extrapolated continuum that dives
    or soars turns the spectrum's edges into spurious lines.
    """
    kx = np.asarray(kx, float)
    return pchip_mono(kx, np.asarray(ky, float))(np.clip(np.asarray(x, float),
                                                         kx.min(), kx.max()))


def sanitise(cont):
    """Blank out zero and non-finite continuum points, which cannot be divided by."""
    cont = np.asarray(cont, float)
    return np.where(np.isfinite(cont) & (cont != 0.0), cont, np.nan)


def local_linear(x, y, percentile=SEED_PERCENTILE, at=None):
    """The straight local pseudo-continuum fitted to ``(x, y)``, evaluated on ``at``
    (default ``x``).

    Uses the same upper-envelope estimator as :func:`seed_local`. Measurement workflows
    keep their nodes fixed; this function is the separate automatic array operation.

    What it is for: the far wings of very strong lines tens of Angstrom outside a narrow
    window depress it as a smooth pseudo-continuum rather than a feature. That depression
    is what a locally renormalised equivalent width divides out, which is also what makes
    a synthesis over a short line-list skirt usable: the lines it omits contribute almost
    only the slope this removes.

    Falls back to a flat :func:`flat_level` when a line cannot be fitted.
    """
    x = np.asarray(x, float)
    at = x if at is None else np.asarray(at, float)
    fit = upper_envelope_fit(x, y, percentile)
    if fit is None:
        level = flat_level(y, percentile) if x.size else 1.0
        return np.full(at.shape, level if np.isfinite(level) and level != 0 else 1.0)
    slope, intercept = fit
    return slope * at + intercept


def divide_local(spec, percentile=SEED_PERCENTILE, window=None):
    """``spec`` divided by its own local straight pseudo-continuum (see :func:`local_linear`).

    The continuum is fitted over ``window`` (in ``spec``'s x unit; its whole extent when
    omitted) and divided out everywhere. For a line that window is the synth window, wider
    than the one the width is integrated over, exactly as the markers are placed; a
    synthesis padded beyond it keeps its padding for convolution.
    """
    x, y = (np.asarray(v, float) for v in spec.values[:2])
    inside = slice(None) if window is None else (x >= window[0]) & (x <= window[1])
    cont = local_linear(x[inside], y[inside], percentile, at=x)
    return spec / np.where(cont == 0, 1.0, cont)


def reject_outliers(spec, *, resolution=None, window=None, percentile=None):
    """The spectrum a continuum spline should be fitted to.

    Low-pass ``spec`` to ``resolution`` (``R = lambda/dlambda``; skipped when unset), then
    drop points whose excess over a running mean taken across ``window`` (in ``spec``'s
    x-unit) exceeds the ``percentile``-th percentile of those excesses. What survives is the
    part of the spectrum that behaves like a continuum; absorption cores and noise spikes do
    not. Fewer than two survivors leaves the spectrum unfiltered rather than empty.
    """
    from scipy.ndimage import uniform_filter1d

    base = spec.lowpass(resolution) if resolution and resolution > 0 else spec
    if not (window and window > 0 and percentile and 0 < percentile < 100):
        return base
    xx, yy = base.values
    dx = abs(float(np.mean(np.diff(xx)))) if len(xx) > 1 else 0.0
    nstep = max(1, int(round(window / dx))) if dx else 1
    excess = yy - uniform_filter1d(yy, nstep, mode="nearest")  # vs the local mean
    mask = excess <= np.percentile(excess, percentile)         # keep points within it
    if mask.sum() < 2:
        return base
    return Spectrum(base.xdata[mask], yy[mask], yunit=base.yunit, xref=base.xref)


def model_continuum(continuum, x, window=None):
    """Evaluate polynomial coefficients or (wavelength, flux) nodes on ``x``.

    Polynomial coefficients follow NumPy order on a dimensionless coordinate spanning
    [-1, 1] across ``window``. Nodes use the same x units as ``x`` and hold their endpoint
    levels outside their span. The result multiplies a normalised model spectrum.
    """
    values = np.asarray(continuum, dtype=float)
    xx = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(values)) or values.size == 0:
        raise ValueError("Continuum values must be finite and nonempty")
    if values.ndim == 2 and values.shape[1] == 2:
        if len(values) < 2 or np.any(np.diff(values[:, 0]) <= 0):
            raise ValueError("Continuum nodes need at least two increasing wavelengths")
        return evaluate(values[:, 0], values[:, 1], xx)
    if values.ndim != 1:
        raise ValueError("Continuum must be a coefficient list or an (N, 2) node array")
    lo, hi = window if window is not None else (xx.min(), xx.max())
    if not hi > lo:
        raise ValueError("A polynomial continuum needs a nonzero wavelength span")
    return np.polyval(values, 2 * (xx - lo) / (hi - lo) - 1)


def seed_local(spec, window):
    """Two continuum nodes at the synth-window edges, seeded from the upper 10% fit."""
    import astropy.units as u
    sub = spec.convert_units(xunit=u.AA).crop(window=window)
    xx, yy = sub.values
    fit = upper_envelope_fit(xx, yy)
    edges = np.asarray(window, dtype=float)
    levels = (np.polyval(fit, edges) if fit is not None
              else np.full(2, flat_level(yy, SEED_PERCENTILE)))
    if not np.all(np.isfinite(levels)) or np.any(levels <= 0):
        raise ValueError("Local continuum levels must be positive and finite")
    return tuple(zip(edges.tolist(), levels.tolist()))
