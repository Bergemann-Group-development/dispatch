"""Fe excitation–ionization balance numerics for the stellar-parameter wizard (pure).

The classic spectroscopic determination of (Teff, logg, [Fe/H], vmic) from per-line Fe
abundances: Teff zeroes the slope of A(Fe I) vs excitation energy, vmic the slope of
A(Fe I) vs reduced EW, logg the Fe I/Fe II ionization offset, and [Fe/H] the consistency
of the mean A(Fe) with the model metallicity. Everything here is numpy-only — the
synthesis-driven orchestration (curve-of-growth inversion per line at the trial
parameters) lives in :mod:`sparamtool`, mirroring the :mod:`abundance` split.

The solver is a damped quasi-Newton fixed point on the four decoupled (diagnostic,
parameter) pairs: iteration 1 uses conservative default sensitivities, later iterations
refine each diagonal Jacobian entry by secant from the observed diagnostic change
(sign-guarded and magnitude-clamped, since the diagonal secant is contaminated by
cross-coupling). Steps are clamped per iteration and to the atmosphere grid's hull.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields

import numpy as np

from abundance import LineResult
from tumult.abundance import reduced_ew

PARAM_KEYS = ("teff", "logg", "feh", "vmic")
DIAG_KEYS = ("s_exc", "s_rew", "d_ion", "d_feh")

# which diagnostic drives which parameter (the classic decoupling)
DIAG_FOR = {"teff": "s_exc", "vmic": "s_rew", "logg": "d_ion", "feh": "d_feh"}

# conservative default sensitivities d(diagnostic)/d(parameter) for an FGK star. The
# signs encode the classic hand rules through the Newton step -d/J: a positive
# excitation slope means Teff is too low (raise it → J_teff < 0), a positive REW slope
# means vmic is too low (raise it → J_vmic < 0), A(Fe II) > A(Fe I) means logg is too
# high (higher logg weakens Fe II lines via the H- continuous opacity, inflating the
# inverted A(Fe II) → J_logg > 0), and the model [Fe/H] cancels out of the consistency
# term almost one-to-one (J_feh ≈ -1).
DEFAULT_JAC = {"teff": -5.0e-5,  # dex/eV per K
               "logg": 0.40,     # dex per dex
               "feh": -0.90,     # dex per dex
               "vmic": -0.15}    # dex/dex(REW) per km/s

# per-iteration step clamps: one bad diagnostic must not fling the solver off the grid
STEP_CLAMP = {"teff": 250.0, "logg": 0.4, "feh": 0.3, "vmic": 0.4}

# secant Jacobian guard: keep the default's sign, clamp the magnitude to this window
# around the default (the diagonal secant picks up cross-terms and noise)
JAC_CLAMP = (0.2, 5.0)

# convergence tolerances on |diagnostic|
DEFAULT_TOLS = {"s_exc": 0.002, "s_rew": 0.01, "d_ion": 0.01, "d_feh": 0.01}


def is_ion_stage(ion, stage):
    """True when a species label like ``"Fe I"`` names ionization ``stage`` (1 → I, 2 → II)."""
    parts = str(ion).split()
    numeral = {1: "I", 2: "II", 3: "III"}.get(stage)
    return len(parts) >= 2 and parts[1] == numeral


def weighted_slope(x, y, sigma=None):
    """Inverse-variance weighted least-squares line fit.

    Returns ``(slope, slope_err, intercept)`` (NaNs when fewer than 3 finite points or
    the x values are degenerate). ``sigma`` (optional) are 1σ errors on ``y``; missing /
    non-finite / non-positive entries fall back to the median finite sigma, or equal
    weights when none is usable.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nan = float("nan")
    if sigma is None:
        s = np.ones_like(y)
    else:
        s = np.asarray(sigma, dtype=float).copy()
        bad = ~np.isfinite(s) | (s <= 0)
        good = s[~bad]
        s[bad] = float(np.median(good)) if good.size else 1.0
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(s)
    x, y, s = x[keep], y[keep], s[keep]
    if x.size < 3 or np.ptp(x) == 0:
        return nan, nan, nan
    w = 1.0 / s**2
    W = w.sum()
    xm = (w * x).sum() / W
    ym = (w * y).sum() / W
    sxx = (w * (x - xm) ** 2).sum()
    slope = float((w * (x - xm) * (y - ym)).sum() / sxx)
    slope_err = float(np.sqrt(1.0 / sxx))
    return slope, slope_err, float(ym - slope * xm)


def _wmean(a, sigma):
    """Inverse-variance weighted mean ``(mean, err)`` with the same sigma fallback as
    :func:`weighted_slope` (equal weights when no line carries a usable error)."""
    a = np.asarray(a, dtype=float)
    if a.size == 0:
        return float("nan"), float("nan")
    s = np.asarray(sigma, dtype=float).copy()
    bad = ~np.isfinite(s) | (s <= 0)
    good = s[~bad]
    s[bad] = float(np.median(good)) if good.size else 1.0
    w = 1.0 / s**2
    W = w.sum()
    mean = float((w * a).sum() / W)
    if a.size > 1:
        # scatter-based error of the mean — robust when the formal sigmas are optimistic
        scatter = np.sqrt((w * (a - mean) ** 2).sum() / W)
        err = float(max(np.sqrt(1.0 / W), scatter / np.sqrt(a.size)))
    else:
        err = float(np.sqrt(1.0 / W))
    return mean, err


def diagnostics(lines, feh, solar_fe):
    """The four balance diagnostics from per-line absolute abundances ``ln.a`` = A(Fe).

    Uses the included lines with finite ``a``. Returns a dict with ``s_exc`` /
    ``s_rew`` (Fe I slopes vs E_low and reduced EW), ``d_ion`` (⟨A(Fe II)⟩ − ⟨A(Fe I)⟩),
    ``d_feh`` (⟨A(Fe)⟩ − (A_solar(Fe) + [Fe/H]_model)), an ``_err`` for each, the means
    ``a_fe1``/``a_fe2``/``a_mean`` and the line counts ``n_fe1``/``n_fe2``. Diagnostics
    without enough lines come back NaN (the solver then holds that parameter).
    """
    inc = [ln for ln in lines if ln.include and np.isfinite(ln.a)]
    fe1 = [ln for ln in inc if is_ion_stage(ln.ion, 1)]
    fe2 = [ln for ln in inc if is_ion_stage(ln.ion, 2)]
    sig = lambda lns: [ln.sigma_a for ln in lns]

    s_exc, s_exc_err, _ = weighted_slope([ln.e_low for ln in fe1],
                                         [ln.a for ln in fe1], sig(fe1))
    s_rew, s_rew_err, _ = weighted_slope([reduced_ew(ln.ew, ln.lam0) for ln in fe1],
                                         [ln.a for ln in fe1], sig(fe1))
    a1, a1_err = _wmean([ln.a for ln in fe1], sig(fe1))
    a2, a2_err = _wmean([ln.a for ln in fe2], sig(fe2))
    am, am_err = _wmean([ln.a for ln in inc], sig(inc))
    d_ion = a2 - a1
    d_ion_err = float(np.hypot(a1_err, a2_err))
    d_feh = am - (float(solar_fe) + float(feh))
    return {
        "s_exc": s_exc, "s_exc_err": s_exc_err,
        "s_rew": s_rew, "s_rew_err": s_rew_err,
        "d_ion": float(d_ion), "d_ion_err": d_ion_err,
        "d_feh": float(d_feh), "d_feh_err": am_err,
        "a_fe1": float(a1), "a_fe2": float(a2), "a_mean": float(am),
        "n_fe1": len(fe1), "n_fe2": len(fe2),
    }


def converged(diags, tols=None):
    """True when every finite diagnostic is inside its tolerance (an all-NaN set is not
    converged — nothing was measured)."""
    tols = {**DEFAULT_TOLS, **(tols or {})}
    vals = [diags.get(k, float("nan")) for k in DIAG_KEYS]
    finite = [v for v in vals if np.isfinite(v)]
    if not finite:
        return False
    return all(abs(diags[k]) <= tols[k] for k in DIAG_KEYS if np.isfinite(diags.get(k, float("nan"))))


@dataclass
class SolverState:
    """The quasi-Newton state carried between iterations.

    ``params`` maps the four axes to their current values; ``jac`` holds the running
    diagonal Jacobian estimates d(diag)/d(param); ``bounds`` (optional) is
    ``{key: (lo, hi)}`` from the atmosphere grid hull.
    """
    params: dict
    jac: dict = field(default_factory=lambda: dict(DEFAULT_JAC))
    bounds: dict = field(default_factory=dict)
    prev_params: dict | None = None
    prev_diags: dict | None = None


def step(state, diags, clamps=None):
    """One damped quasi-Newton update; returns the new params dict (and mutates ``state``).

    Each parameter moves by ``-diag / J`` along its own diagnostic, clamped per
    iteration and to the grid bounds. Before stepping, each Jacobian entry is refined by
    secant from the previous iteration's (param, diagnostic) pair when that parameter
    actually moved — guarded to the default's sign and a magnitude window, because the
    diagonal secant absorbs cross-coupling between the axes. A NaN diagnostic (e.g. no
    Fe II lines survived) holds its parameter.
    """
    clamps = {**STEP_CLAMP, **(clamps or {})}
    if state.prev_params is not None and state.prev_diags is not None:
        for key in PARAM_KEYS:
            dk = DIAG_FOR[key]
            dp = state.params[key] - state.prev_params[key]
            dd = diags.get(dk, float("nan")) - state.prev_diags.get(dk, float("nan"))
            if not (np.isfinite(dp) and np.isfinite(dd)) or abs(dp) < 1e-12:
                continue
            sec = dd / dp
            ref = DEFAULT_JAC[key]
            if not np.isfinite(sec) or sec * ref <= 0:
                continue                      # wrong sign / degenerate — keep the estimate
            lo, hi = JAC_CLAMP
            state.jac[key] = float(np.clip(abs(sec), lo * abs(ref), hi * abs(ref))) * np.sign(ref)

    new = dict(state.params)
    for key in PARAM_KEYS:
        d = diags.get(DIAG_FOR[key], float("nan"))
        j = state.jac.get(key, DEFAULT_JAC[key])
        if not (np.isfinite(d) and np.isfinite(j)) or j == 0:
            continue
        delta = float(np.clip(-d / j, -clamps[key], clamps[key]))
        val = new[key] + delta
        if key in state.bounds:
            lo, hi = state.bounds[key]
            val = float(np.clip(val, lo, hi))
        new[key] = val

    state.prev_params = dict(state.params)
    state.prev_diags = dict(diags)
    state.params = new
    return new


def uncertainties(diags, jac):
    """Parameter 1σ from the converged diagnostics: ``σ_p = σ_diag / |J_pp|``.

    First-order and diagonal-only — the honest formal error of the balance method
    itself, not a full covariance. NaN where the diagnostic error (or Jacobian) is
    unusable.
    """
    out = {}
    for key in PARAM_KEYS:
        dk = DIAG_FOR[key]
        err = diags.get(f"{dk}_err", float("nan"))
        j = jac.get(key, float("nan"))
        out[key] = float(abs(err / j)) if np.isfinite(err) and np.isfinite(j) and j != 0 \
            else float("nan")
    return out


@dataclass
class SParamResult:
    """One stellar-parameter determination: the curated Fe line set plus the solve record.

    ``lines`` mixes Fe I and Fe II (distinguished by ``ln.ion``); ``params``/``params0``
    are the solved / initial ``{teff, logg, feh, vmic}``; ``history`` records one dict
    per iteration (``{"params": …, "diags": …}``); ``sigmas`` the propagated parameter
    errors; ``diags`` the final diagnostics. ``stage`` walks "lines" → "solved".
    """
    lines: list                        # list[LineResult]
    params: dict = field(default_factory=dict)
    params0: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    sigmas: dict = field(default_factory=dict)
    diags: dict = field(default_factory=dict)
    jac: dict = field(default_factory=dict)   # final Jacobian diag, for re-aggregation
    converged: bool = False
    solar_fe: float = float("nan")     # A_solar(Fe) used by d_feh (kept for re-aggregation)
    stage: str = "lines"
    profile: str = "gauss"
    snr: float | None = None
    cont_pct: float = float("nan")
    uid: str = ""
    name: str = ""

    def to_dict(self):
        d = {}
        for f in fields(self):
            d[f.name] = [ln.to_dict() for ln in self.lines] if f.name == "lines" \
                else getattr(self, f.name)
        if hasattr(self, "ref_uid"):
            d["ref_uid"] = self.ref_uid
        return d

    @classmethod
    def from_dict(cls, d):
        d = dict(d or {})
        names = {f.name for f in fields(cls)}
        lines = [LineResult.from_dict(x) for x in (d.get("lines") or [])]
        res = cls(lines=lines, **{k: v for k, v in d.items() if k in names and k != "lines"})
        if "ref_uid" in d:
            res.ref_uid = d["ref_uid"]
        return res
