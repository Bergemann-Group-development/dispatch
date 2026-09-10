"""Internal numerical building blocks for :mod:`tumult.abundance`.

The clearly-general pieces of the abundance pipeline live here so they can be reused outside the
TumultStudio2 GUI: curve-of-growth inversion, the inverse-variance element mean, the
fitted-profile FWHM, abundance-scale conversions (A(X) <-> [X/Fe]), the solar reference table
reader, and the EW-error -> abundance-error propagation. The numerics here are numpy-only — no
synthesis engine or GUI dependencies. The reduced equivalent width now lives on
:mod:`tumult.profile` (alongside :attr:`Profile.rew`) and is re-exported here for back-compat.
The synthesis-driven orchestration (screening, per-line fits, gradients) stays in the
application layer, which re-exports these for back-compat.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .const import CC_KMS
from .profile import get_profile_class, reduced_ew  # noqa: F401 (reduced_ew re-export)
from .utils import pchip_mono


def cog_invert(a_grid, ew_grid, ew_meas, ew_err=0.0):
    """Interpolate an abundance from a measured EW on a curve of growth.

    Parameters
    ----------
    a_grid, ew_grid : array_like
        The curve of growth: abundances and the synthetic EW each produces (same units
        for ``ew_grid`` and ``ew_meas``). Need not be sorted; monotonic EW(A) assumed.
    ew_meas : float
        The measured equivalent width.
    ew_err : float, optional
        1σ uncertainty on ``ew_meas``. Propagated to ``sigma_a`` via the local slope.

    Returns
    -------
    dict with keys
        ``a`` : interpolated abundance at ``ew_meas``.
        ``sigma_a`` : ``ew_err * |dA/dEW|`` at the measured point (∞-prone on saturation).
        ``dadew`` : the local ``dA/dEW`` slope (1/EW-unit).
        ``out_of_range`` : True if ``ew_meas`` falls outside the synthesised CoG span.
    """
    a = np.asarray(a_grid, dtype=float)
    w = np.asarray(ew_grid, dtype=float)
    if a.size < 2 or a.size != w.size:
        raise ValueError("cog_invert needs matching a_grid/ew_grid of length >= 2.")
    order = np.argsort(w)
    w, a = w[order], a[order]

    # monotone piecewise-cubic (PCHIP) interpolant A(EW): shape-preserving so it neither
    # overshoots nor introduces spurious extrema on a monotonic CoG, and its analytic
    # derivative gives a smooth local dA/dEW slope for the error propagation.
    spline = pchip_mono(w, a)
    a_val = float(spline(ew_meas))
    dadew = float(spline.derivative()(ew_meas))
    sigma_a = float(abs(ew_err * dadew)) if np.isfinite(dadew) else np.inf
    return {
        "a": a_val,
        "sigma_a": sigma_a,
        "dadew": float(dadew),
        "out_of_range": bool(ew_meas < w[0] or ew_meas > w[-1]),
    }


def aggregate(lines):
    """Inverse-variance weighted element mean over the included lines, with a budget split.

    ``lines`` is any iterable of objects exposing ``include`` (bool), ``a``, ``sigma_a``,
    ``sigma_a_phot`` and ``sigma_a_cont``. Returns a dict: ``a_mean``, ``sigma_mean`` (formal
    inverse-variance error), ``sigma_scatter`` (weighted line-to-line scatter, a robustness
    cross-check), ``sigma_phot`` / ``sigma_cont`` (photon vs continuum shares, summing in
    quadrature to ``sigma_mean``), and ``n`` (lines used). Empty/zero-weight input → NaNs, ``n=0``.
    """
    inc = [ln for ln in lines if ln.include and np.isfinite(ln.sigma_a) and ln.sigma_a > 0]
    if not inc:
        nan = float("nan")
        return {"a_mean": nan, "sigma_mean": nan, "sigma_scatter": nan,
                "sigma_phot": nan, "sigma_cont": nan, "n": 0}
    a = np.array([ln.a for ln in inc], dtype=float)
    s = np.array([ln.sigma_a for ln in inc], dtype=float)
    w = 1.0 / s**2
    wsum = w.sum()
    a_mean = float((w * a).sum() / wsum)
    sigma_mean = float(np.sqrt(1.0 / wsum))
    # weighted line-to-line scatter (population, not corrected — a conservative cross-check)
    sigma_scatter = float(np.sqrt((w * (a - a_mean) ** 2).sum() / wsum))

    # decompose: with weights w_i = 1/σ_i², the mean's variance from a component k is
    # Σ w_i² σ_{k,i}² / (Σ w_i)². Photon and continuum shares add to sigma_mean² exactly
    # when σ_i² = σ_phot,i² + σ_cont,i².
    sphot = np.array([ln.sigma_a_phot for ln in inc], dtype=float)
    scont = np.array([ln.sigma_a_cont for ln in inc], dtype=float)
    denom = wsum**2
    sigma_phot = float(np.sqrt((w**2 * sphot**2).sum() / denom))
    sigma_cont = float(np.sqrt((w**2 * scont**2).sum() / denom))
    return {"a_mean": a_mean, "sigma_mean": sigma_mean, "sigma_scatter": sigma_scatter,
            "sigma_phot": sigma_phot, "sigma_cont": sigma_cont, "n": len(inc)}


def aggregate_unweighted(lines):
    """Unweighted mean + population line-to-line scatter of the included values.

    The companion to :func:`aggregate` for measurements that carry **no usable per-line
    formal error**, where inverse-variance weighting is undefined and the line-to-line
    dispersion is the only uncertainty estimate. Note that *fitted* profiles always
    carry a covariance-derived EW error (residual-scaled, i.e. model mismatch — nonzero
    even on noiseless synthetic spectra), so results measured through them can and
    should use the weighted :func:`aggregate`; this fallback remains for trapezoidal /
    pchip integrations without a noise model and for results restored from sessions
    predating the error propagation. ``lines`` expose ``include`` and ``a``.
    Returns ``{a_mean, sigma_scatter, n}`` (NaNs / ``n=0`` when nothing is included).
    """
    a = np.array([ln.a for ln in lines
                  if ln.include and np.isfinite(ln.a)], dtype=float)
    if a.size == 0:
        return {"a_mean": float("nan"), "sigma_scatter": float("nan"), "n": 0}
    return {"a_mean": float(a.mean()),
            "sigma_scatter": float(a.std(ddof=1)) if a.size > 1 else 0.0,
            "n": int(a.size)}


# --- abundance scales (A(X) <-> [X/Fe]) and the solar reference table --------------------


def xfe_to_abund(xfe, el, feh, solar):
    """Convert ``[X/Fe]`` to an absolute ``A(X)``: ``A(X) = A_solar(X) + [Fe/H] + [X/Fe]``.

    Exact for metals with no alpha/r-process enhancement (the configuration used here):
    MULTI3D scales every element except H/He/Li by [Fe/H], so a relative override equals
    [X/Fe]. Returns ``xfe`` unchanged if the element is missing from ``solar``.
    """
    a0 = solar.get(el)
    return xfe if a0 is None else xfe + a0 + feh


def abund_to_xfe(a, el, feh, solar):
    """Inverse of :func:`xfe_to_abund`: ``[X/Fe] = A(X) - A_solar(X) - [Fe/H]``."""
    a0 = solar.get(el)
    return a if a0 is None else a - a0 - feh


def convert_abundances(text, to_relative, feh, solar):
    """Convert an ``El=val ...`` override string between absolute A(X) and [X/Fe].

    ``to_relative=True``: the values are currently absolute A(X); emit ``[X/Fe] =
    A(X) - A_solar(X) - [Fe/H]``. ``to_relative=False``: invert. Whitespace- or
    comma-separated; terms whose element is unknown (or unparsable) pass through unchanged.
    """
    out = []
    for term in re.split(r"[\s,]+", str(text).strip()):
        if not term:
            continue
        el, sep, val = term.partition("=")
        a0 = solar.get(el)
        if not sep or a0 is None:
            out.append(term)
            continue
        try:
            x = float(val)
        except ValueError:
            out.append(term)
            continue
        x = x - a0 - feh if to_relative else x + a0 + feh
        out.append(f"{el}={x:g}")
    return " ".join(out)


def read_solar_abundances(path):
    """``{element: A(X)}`` parsed from a solar-reference abundance file.

    The format is the dispatch ``abund.txt`` convention: one ``SYMBOL VALUE`` pair per
    line (extra columns ignored, unparsable lines skipped, symbols capitalised).
    Returns ``{}`` if the file is unreadable.
    """
    table = {}
    try:
        lines = Path(path).expanduser().read_text().splitlines()
    except OSError:
        return {}
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0].isalpha():
            try:
                table[parts[0].capitalize()] = float(parts[1])
            except ValueError:
                continue
    return table


# --- EW-error -> abundance-error propagation ----------------------------------------------


def velocity_hwhm(fwhm_aa, lam0_aa):
    """Half width at half maximum expressed as a velocity [km/s] (``HWHM_λ / λ0 · c``).

    Lets a vshift screen compare a fitted radial-velocity shift (km/s) against the line's
    own width. Non-finite / non-positive inputs give ``NaN`` so the screen is skipped.
    """
    if not (np.isfinite(fwhm_aa) and np.isfinite(lam0_aa)) or lam0_aa <= 0:
        return float("nan")
    return 0.5 * abs(fwhm_aa) / lam0_aa * CC_KMS


def sigma_a_from_rew(prof, drewda):
    """Convert a measured line's EW errors to an abundance error through the gradient ``dREW/dA``.

    Works in reduced-EW space (unit independent): the :class:`tumult.profile.Profile` already
    carries ``σ_REW = σ_W / (W ln10)`` per component (``rew_err_phot`` / ``rew_err_cont``), and
    ``σ_A = σ_REW / |dREW/dA|``. Returns ``(sigma_a, sigma_a_phot, sigma_a_cont)``; a non-finite or
    zero gradient / EW blows the error up to ``inf`` (as a saturated or vanished line should).
    """
    g = abs(float(drewda))
    if not np.isfinite(prof.maa) or prof.maa <= 0 or not np.isfinite(g) or g == 0:
        return np.inf, np.inf, np.inf
    s_phot = prof.rew_err_phot / g
    s_cont = prof.rew_err_cont / g
    return float(np.hypot(s_phot, s_cont)), float(s_phot), float(s_cont)


def profile_fwhm(pout, func, profile):
    """Full width at half depth of a fitted line profile, in the fit's wavelength unit.

    Thin shim over :attr:`tumult.profile.FittedProfile.fwhm` (Gaussian widths analytic, the
    Voigt measured numerically off the fitted ``func``). ``pout`` are the fit parameters and
    ``func`` the fitted ``f(λ)``.
    """
    return get_profile_class(profile)(None, params=pout, func=func).fwhm
