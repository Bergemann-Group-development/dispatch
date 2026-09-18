"""Synthetic spectra whose properties a reader can work out by hand.

Everything here is an *input* factory. Nothing calls the code under test, so a case
built on one of these may compare against the closed form of the shape below.

Conventions shared by every factory:

* wavelengths in Angstrom, on the air scale that the user-facing side of OpenMULTI3D
  uses (see *wavelength medium* in ``docs/GLOSSARY.md``). No conversion is applied
  here, and a measurement of width or centre is unaffected by the choice, since both
  are read off the same axis they were built on;
* flux already continuum normalised, so the continuum is exactly 1.0 everywhere and
  no continuum placement enters a measurement made on it;
* a uniform wavelength step, so trapezoidal weights are constant and the sampling can
  be quoted in samples per sigma.
"""
from __future__ import annotations

import numpy as np

# The reference line: an absorption feature deep enough to measure and narrow enough
# that a 4 A window holds many widths of it.
LINE_DEPTH = 0.42        # fractional depth at the centre, on a unit continuum
LINE_SIGMA = 0.11        # Gaussian standard deviation [AA]
LINE_CENTRE = 5000.05    # true centre [AA]
LINE_LAM0 = 5000.0       # the catalogue wavelength a measurement is asked for [AA]


def gaussian_line(n_points=801, half_span=2.0, *, depth=LINE_DEPTH, sigma=LINE_SIGMA,
                  centre=LINE_CENTRE, lam0=LINE_LAM0):
    """``(x, y)`` for a single Gaussian absorption line on a unit continuum.

    ``y = 1 - depth * exp(-0.5 ((x - centre)/sigma)^2)`` over
    ``[lam0 - half_span, lam0 + half_span]`` with ``n_points`` uniform samples.

    The defaults give a step of ``4/800 = 0.005 AA``, that is 22 samples per sigma, and
    a half window of ``2/0.11 = 18 sigma``: the profile is below ``1e-70`` of its depth
    at the edges, so the window truncates nothing measurable.

    The line is one unblended feature in the weak-line regime, with no saturation and no
    second component, which is what makes its equivalent width the closed form
    ``W = depth * sqrt(2 pi) * sigma``.
    """
    x = np.linspace(lam0 - half_span, lam0 + half_span, n_points)
    y = 1.0 - depth * np.exp(-0.5 * ((x - centre) / sigma) ** 2)
    return x, y


def gaussian_line_equivalent_width_maa(depth=LINE_DEPTH, sigma=LINE_SIGMA):
    """Closed-form equivalent width of :func:`gaussian_line`, in milli-Angstrom.

    ``W = integral (1 - y) dx = depth * sqrt(2 pi) * sigma`` in Angstrom, times 1000.
    This is the analytic value; it is never obtained by measuring the line.
    """
    return depth * np.sqrt(2.0 * np.pi) * sigma * 1e3


def flux_centroid(x, y):
    """Absorption-weighted centroid of ``(x, y)`` on a unit continuum [same unit as x].

    ``sum((1 - y) x) / sum(1 - y)``. On a uniform grid this is the centre of a symmetric
    profile exactly, which is what makes it an independent estimate of a shifted line's
    position.
    """
    weight = 1.0 - np.asarray(y, dtype=float)
    return float((weight * np.asarray(x, dtype=float)).sum() / weight.sum())
