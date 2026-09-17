"""Excitation–ionization balance numerics (pure).

The classic spectroscopic determination of (Teff, logg, [Fe/H], vmic) from per-line
abundances of one element in two ionization stages: Teff zeroes the slope of A(X I) vs
excitation energy, vmic the slope of A(X I) vs reduced EW, logg the X I/X II ionization
offset, and the metallicity the consistency of the mean A(X) with the model's. Everything
in the first half is numpy-only; :func:`ionbalance` at the end drives the syntheses that
turn measured EWs into the per-line abundances those numerics consume.

Nothing here is Fe-specific. For an element other than Fe the metallicity axis closes only
under a scaled-solar
assumption — ``d_feh`` compares ⟨A(X)⟩ with ``A_solar(X) + [Fe/H]_model``, so the caller is
asserting [X/Fe] = 0. That is exactly right for Fe, defensible for most iron-peak elements
and known to be wrong for the α elements at low metallicity; the choice belongs to the
caller, which is why the solar abundance arrives as an argument.

:func:`solve_coupled` steps the free parameters on the free diagnostics as one residual
vector with a full Jacobian, so a Teff step that moves the ionisation offset is accounted
for. :func:`ionbalance` computes that Jacobian from the same node runs as the diagnostics
(the slope of each trial star's interpolation cell, see :func:`ionbalance`). Steps are clamped per iteration and to the atmosphere grid's
hull.

A parameter may be held *fixed*: it then never steps, and its diagnostic is still computed
and reported but does not gate convergence — the residual at a fixed parameter is the
interesting output ("at your log g the ionization balance is off by 0.08 dex"). Fixing all
four turns the solve into a single diagnostic evaluation.
"""
from __future__ import annotations

import logging
from pathlib import Path
from dataclasses import dataclass

import numpy as np

from .abundance import reduced_ew, xfe_to_abund

PARAM_KEYS = ("teff", "logg", "feh", "vmic")
DIAG_KEYS = ("s_exc", "s_rew", "d_ion", "d_feh")

# which diagnostic drives which parameter (the classic decoupling)
DIAG_FOR = {"teff": "s_exc", "vmic": "s_rew", "logg": "d_ion", "feh": "d_feh"}
# sign of d DIAG_FOR[p] / d p: hotter flattens the excitation slope, more vmic
# desaturates strong lines, higher log g raises Fe II against Fe I, and a model
# [Fe/H] step moves the mean abundance by less than itself
DIAG_SIGN = {"teff": -1, "vmic": -1, "logg": 1, "feh": -1}

# per-iteration step clamps: one bad diagnostic must not fling the solver off the grid
STEP_CLAMP = {"teff": 250.0, "logg": 0.4, "feh": 0.3, "vmic": 0.4}

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


def diagnostics(lines, feh, solar_x):
    """The four balance diagnostics from per-line absolute abundances ``ln.a`` = A(X).

    Uses the included lines with finite ``a``. Returns a dict with ``s_exc`` /
    ``s_rew`` (neutral-stage slopes vs E_low and reduced EW), ``d_ion``
    (⟨A(X II)⟩ − ⟨A(X I)⟩), ``d_feh`` (⟨A(X)⟩ − (``solar_x`` + [Fe/H]_model), i.e. the
    scaled-solar metallicity consistency), an ``_err`` for each, the means
    ``a_1``/``a_2``/``a_mean`` and the line counts ``n_1``/``n_2``. Diagnostics
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
    d_feh = am - (float(solar_x) + float(feh))
    return {
        "s_exc": s_exc, "s_exc_err": s_exc_err,
        "s_rew": s_rew, "s_rew_err": s_rew_err,
        "d_ion": float(d_ion), "d_ion_err": d_ion_err,
        "d_feh": float(d_feh), "d_feh_err": am_err,
        "a_1": float(a1), "a_2": float(a2), "a_mean": float(am),
        "n_1": len(fe1), "n_2": len(fe2),
    }


def free_diags(free=None):
    """The diagnostic keys that gate convergence: those driving a free parameter.

    ``free`` is ``{param: bool}``; ``None`` means every parameter is free. A fixed
    parameter's diagnostic is still computed and reported — it just cannot be nulled, so
    letting it gate would make such a solve unconvergeable by construction.
    """
    if free is None:
        return list(DIAG_KEYS)
    return [DIAG_FOR[k] for k in PARAM_KEYS if free.get(k, True)]


def converged(diags, tols=None, free=None):
    """True when every free diagnostic is finite and inside its tolerance."""
    tols = {**DEFAULT_TOLS, **(tols or {})}
    keys = free_diags(free)
    return bool(keys) and all(
        np.isfinite(diags.get(k, float("nan"))) and abs(diags[k]) <= tols[k]
        for k in keys)


def balance_status(diags, tols=None, free=None):
    """Classify an evaluated balance independently of its iteration budget."""
    keys = free_diags(free)
    if not keys:
        return "diagnostic"
    if any(not np.isfinite(diags.get(k, float("nan"))) for k in keys):
        return "undetermined"
    return "converged" if converged(diags, tols, free) else "iteration_limit"


# --- the synthesis-driven solve ------------------------------------------------------------
# Everything above is numpy on measured lines. What follows drives a synthesis: it inverts
# each line's EW to an abundance through a curve of growth computed at the trial parameters,
# turns those abundances into the four diagnostics above, and steps. It is here rather than
# in the GUI because nothing in it is about rendering — a script, a notebook or an MCP
# adapter runs the same solve by handing it a request and a measurement.

BLEND_MAX = 0.30           # blend screen: contaminant EW fraction ceiling
VMIC_BOUNDS = (0.1, 5.0)   # vmic is a synthesis input, not a grid axis — keep it sane


def entry_star(params):
    """The reference atmosphere a line identification ranks strengths on."""
    from . import refstars
    return refstars.nearest(params["teff"], params["logg"], params["feh"])


def identify_catalog(star, lo, hi, linelist, atmosphere=None):
    """``(frame, note)``: the transitions in ``[lo, hi]`` Å ranked for ``star``.

    HOP strengths on a model atmosphere are what makes the identification a statement about
    *this* star: on ``atmosphere`` (a MARCS model) when given, else on the reference
    atmosphere ``star`` resolves to through the registered MARCS grid. The master line list
    is the fallback when HOP cannot be run, and then the note says so, because a log
    gf-ordered identification is a much weaker claim.
    """
    from . import linestrength
    from .linelist import load_master
    try:
        if atmosphere is not None:
            frame = linestrength.band_on_atmosphere(atmosphere, str(linelist), lo, hi)
            label = f"HOP strengths on {Path(atmosphere).name}"
        else:
            frame = linestrength.compute_band(star, str(linelist),
                                              *linestrength.band_bounds(lo, hi))
            label = f"HOP strengths on {star.label}"
        if frame is not None and len(frame):
            return frame, label
    except Exception:
        logging.exception("Line-strength band failed; falling back to the line list")
    try:
        return load_master(linelist), "line-list order (no HOP strengths)"
    except Exception:
        logging.exception("Master line list unreadable")
        return None, "no line list — identities unknown"


def measure_balance_lines(result, masks, spec, linelist, atmosphere=None):
    """Fill ``result.lines`` from the windows, identified on the entry point's atmosphere.

    ``atmosphere`` names that MARCS model; without it the reference star nearest the entry
    point is used (see :func:`identify_catalog`).

    Raises :class:`ValueError` when no window could be identified as either stage of the
    element: the balance reads the *stage* off this identification, and the mask set's
    element cannot supply one. Windows that are dominated by another species, or whose
    dominant transition is the third stage or above, are switched off rather than dropped —
    the flag stays visible and the caller's toggle overrides it.
    """
    from . import linediag
    star = entry_star(result.params0) if atmosphere is None else None
    lo = min(m.lmin for m in masks)
    hi = max(m.lmax for m in masks)
    catalog, note = identify_catalog(star, lo, hi, linelist, atmosphere)
    lines = linediag.measure_masks(spec, masks, profile=result.profile, snr=result.snr,
                                   cont_pct=result.cont_pct, catalog=catalog,
                                   element=result.element)
    result.ident_star = star.label if star is not None else Path(atmosphere).name
    result.ident_note = note

    for ln in lines:
        stage = is_ion_stage(ln.ion, 1) or is_ion_stage(ln.ion, 2)
        if "off_element" in ln.flags or not stage:
            ln.include = False
            if not stage and "off_element" not in ln.flags:
                ln.flags = ln.flags + ("no_stage",)
    result.lines = lines
    n1, n2 = result.counts()
    if not n1 and not n2:
        raise ValueError(
            f"No window in “{result.maskset_name}” could be identified as "
            f"{result.element} I or {result.element} II — the balance needs the stage.")
    result.name = result.label()
    return result


def grid_bounds(grid, atmos_file=None):
    """Box bounds of the atmosphere grid, as a hull proxy.

    The interpolation itself raises on a truly unreachable query, which the solve surfaces;
    these bounds only keep a step from walking off the grid in the first place.
    """
    from .synthesis.snapshot_grid import discover_grid
    if atmos_file or not grid.get("folder"):
        return {"vmic": VMIC_BOUNDS}
    try:
        files = discover_grid(grid)
    except Exception:
        return {"vmic": VMIC_BOUNDS}
    keys = np.array([k[:3] for k in files], dtype=float)
    bounds = {name: (float(keys[:, i].min()), float(keys[:, i].max()))
              for i, name in enumerate(("teff", "logg", "feh"))}
    vmics = [k[3] for k in files if len(k) > 3]
    bounds["vmic"] = (float(min(vmics)), float(max(vmics))) if vmics else VMIC_BOUNDS
    return bounds


# --- the coupled solve -----------------------------------------------------------------------
# Numerical contract of solve_coupled:
#
# * residuals: the free diagnostics, each divided by its tolerance, so every component is
#   dimensionless and "converged" means every |r_k| <= 1;
# * variables: each free parameter divided by its finite-difference step (FD_STEPS), so the
#   damping treats one step of every parameter alike;
# * Jacobian: the one ``evaluate`` returns with each point, whenever it is finite;
#   otherwise forward differences at the point, one evaluation per free parameter at
#   :func:`probe_point`, then updated by Broyden's rank-one formula after an accepted step
#   and recomputed where a step is rejected on a Broyden estimate;
# * step: Levenberg-Marquardt, (A^T A + lambda diag(A^T A)) dx = -A^T r, with lambda
#   starting at LM_LAMBDA, divided by 3 on an accepted step and multiplied by 4 on a
#   rejected one (a step is accepted when it lowers |r|); the step is scaled down as a
#   whole to respect STEP_CLAMP and then projected onto the grid bounds;
# * termination: "converged" (every free diagnostic within tolerance), "undetermined"
#   (a free diagnostic or Jacobian column is not finite), "boundary" (the projected step
#   is zero), "iteration_limit" (max_iter trial steps), "diagnostic" (nothing free).
#   Finite-difference evaluations do not count towards max_iter.

FD_STEPS = {"teff": 100.0, "logg": 0.1, "feh": 0.1, "vmic": 0.2}
LM_LAMBDA = 1e-2
COUPLED_MAX_ITER = 12
PROBE_HALVINGS = 4


def probe_point(params, key, step, bounds=None):
    """``params`` moved by ``step`` along ``key`` for a one-sided difference, or None.

    The step goes forward unless that leaves ``bounds[key]``, then backward; it is clipped
    to the bounds, and None is returned when no move remains.
    """
    lo, hi = (bounds or {}).get(key, (-np.inf, np.inf))
    value = params[key] + step
    if value > hi:
        value = params[key] - step
    value = float(np.clip(value, lo, hi))
    if value == params[key]:
        return None
    return {**params, key: value}


@dataclass
class CoupledOutcome:
    """What :func:`solve_coupled` found: the accepted point and how the solve ended."""
    params: dict
    diags: dict
    state: object
    history: list
    termination: str
    jacobian: dict
    sigmas: dict
    covariance: dict

    @property
    def converged(self):
        return self.termination == "converged"


def coupled_uncertainties(diags, jacobian, free=None):
    """Parameter 1σ and covariance from the diagnostic errors through the full Jacobian.

    ``jacobian`` is ``{diagnostic: {parameter: d diagnostic / d parameter}}`` over the free
    pairs. With diagonal diagnostic errors ``C_d`` the parameter covariance is
    ``(J^T C_d^-1 J)^-1``. NaN where an error or the inversion is unusable and for fixed
    parameters.
    """
    keys = [k for k in PARAM_KEYS if free is None or free.get(k, True)]
    nan = float("nan")
    sigmas = {k: nan for k in PARAM_KEYS}
    covariance = {}
    dkeys = [DIAG_FOR[k] for k in keys]
    try:
        J = np.array([[float(jacobian[d][k]) for k in keys] for d in dkeys])
        err = np.array([float(diags.get(f"{d}_err", nan)) for d in dkeys])
        if not keys or not (np.all(np.isfinite(J)) and np.all(np.isfinite(err))
                            and np.all(err > 0)):
            return sigmas, covariance
        cov = np.linalg.inv(J.T @ np.diag(1.0 / err ** 2) @ J)
    except (KeyError, TypeError, ValueError, np.linalg.LinAlgError):
        return sigmas, covariance
    for i, k in enumerate(keys):
        sigmas[k] = float(np.sqrt(cov[i, i])) if cov[i, i] >= 0 else nan
        covariance[k] = {q: float(cov[i, j]) for j, q in enumerate(keys)}
    return sigmas, covariance


def solve_coupled(evaluate, params0, *, free=None, bounds=None, tols=None,
                  max_iter=COUPLED_MAX_ITER, steps=None, clamps=None, report=None,
                  cancel_token=None):
    """Null the free balance diagnostics with a damped, coupled Newton iteration.

    ``evaluate(params) -> (diags, state, jacobian)`` computes the diagnostics at a
    parameter point; ``state`` is whatever the caller wants back for the accepted point
    (the per-line abundances) and ``jacobian`` ``{diagnostic: {parameter: d/dp}}`` at it,
    or None to have it differenced here. ``free`` maps parameters to whether they are fitted, ``bounds`` gives
    ``{param: (lo, hi)}``, ``tols`` overrides :data:`DEFAULT_TOLS`, ``steps``
    :data:`FD_STEPS` and ``clamps`` :data:`STEP_CLAMP`. See the contract above this
    function. ``report(history)`` is called after every evaluation.
    """
    from .synthesis.cancel import CancelledError

    free = {k: bool((free or {}).get(k, True)) for k in PARAM_KEYS}
    keys = [k for k in PARAM_KEYS if free[k]]
    dkeys = [DIAG_FOR[k] for k in keys]
    tol = {**DEFAULT_TOLS, **(tols or {})}
    scale = np.array([float({**FD_STEPS, **(steps or {})}[k]) for k in keys])
    clamp = np.array([float({**STEP_CLAMP, **(clamps or {})}[k]) for k in keys])
    bounds = dict(bounds or {})
    report = report or (lambda _history: None)
    history = []

    def project(params):
        out = dict(params)
        for k in PARAM_KEYS:
            if k in bounds:
                lo, hi = bounds[k]
                out[k] = float(np.clip(out[k], lo, hi))
        return out

    def run(params, kind):
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        diags, state, jacobian = evaluate(dict(params))
        history.append({"params": dict(params), "diags": dict(diags), "kind": kind})
        report(history)
        return diags, state, jacobian

    def scaled(jacobian):
        """``jacobian`` as d r / d x (NaN where missing), or None without one."""
        if jacobian is None:
            return None
        return np.array([[float((jacobian.get(d) or {}).get(k, np.nan)) * scale[j] / tol[d]
                          for j, k in enumerate(keys)] for d in dkeys])

    def residual(diags):
        return np.array([float(diags.get(d, float("nan"))) / tol[d] for d in dkeys])

    def fd_jacobian(params, r):
        """Columns d r / d x at ``params``, x the parameters in units of their step."""
        J = np.empty((len(dkeys), len(keys)))
        for j, k in enumerate(keys):
            probe = probe_point(params, k, scale[j], bounds)
            if probe is None:
                J[:, j] = np.nan
                continue
            dx = (probe[k] - params[k]) / scale[j]
            diags, _state, _jacobian = run(probe, "jacobian")
            J[:, j] = (residual(diags) - r) / dx
        return J

    def outcome(termination, J):
        jacobian = {}
        if J is not None and np.all(np.isfinite(J)):
            # back to diagnostic units per parameter unit
            raw = J * np.array([tol[d] for d in dkeys])[:, None] / scale[None, :]
            jacobian = {d: {k: float(raw[i, j]) for j, k in enumerate(keys)}
                        for i, d in enumerate(dkeys)}
        sigmas, covariance = coupled_uncertainties(diags, jacobian, free) \
            if jacobian else ({k: float("nan") for k in PARAM_KEYS}, {})
        return CoupledOutcome(dict(params), dict(diags), state, history, termination,
                              jacobian, sigmas, covariance)

    params = project({k: float(v) for k, v in params0.items()})
    diags, state, given = run(params, "step")
    if not keys:
        return outcome("diagnostic", None)
    J, fresh, lam = None, False, LM_LAMBDA
    for _it in range(int(max_iter)):
        r = residual(diags)
        if not np.all(np.isfinite(r)):
            return outcome("undetermined", J)
        if converged(diags, tol, free):
            return outcome("converged", scaled(given) if J is None else J)
        if J is None:
            J = scaled(given)
            if J is None or not np.all(np.isfinite(J)):
                J = fd_jacobian(params, r)
            fresh = True
            if not np.all(np.isfinite(J)):
                return outcome("undetermined", J)
        H, g = J.T @ J, J.T @ r
        damped = H + lam * np.diag(np.maximum(np.diag(H), 1e-12))
        dx = -np.linalg.lstsq(damped, g, rcond=None)[0]
        dp = dx * scale
        over = np.max(np.abs(dp) / clamp)
        if over > 1:
            dp = dp / over
        trial = project({**params, **{k: params[k] + dp[j] for j, k in enumerate(keys)}})
        moved = np.array([trial[k] - params[k] for k in keys])
        if not np.any(np.abs(moved) > 1e-12):
            return outcome("boundary", J)
        new_diags, new_state, new_given = run(trial, "step")
        r_new = residual(new_diags)
        if np.all(np.isfinite(r_new)) and np.linalg.norm(r_new) < np.linalg.norm(r):
            if new_given is None:
                dxm = moved / scale
                J = J + np.outer((r_new - r) - J @ dxm, dxm) / float(dxm @ dxm)
            else:
                J = None          # the accepted point brought its own
            params, diags, state, given, fresh = trial, new_diags, new_state, new_given, False
            lam = max(lam / 3.0, 1e-8)
        else:
            history[-1]["kind"] = "rejected"
            lam *= 4.0
            if not fresh:
                J = None
    J = scaled(given) if J is None else J
    r = residual(diags)
    if not np.all(np.isfinite(r)):
        return outcome("undetermined", J)
    return outcome("converged" if converged(diags, tol, free) else "iteration_limit", J)


def ionbalance(request, *, target, line_mask=None, measured=None, method="cog",
               background_lines=True, free=None, tols=None, max_iter=COUPLED_MAX_ITER,
               profile="voigt", snr=None, cont_pct=0.0, limit=1.0, max_rounds=8, steps=None,
               blend_max=BLEND_MAX, weak_frac=0.10, progress=None, cancel_token=None):
    """Infer stellar parameters from the excitation-ionisation balance of ``target``.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest` on an atmosphere grid; its
    star is the solver's starting point and the atmosphere the line stages are identified
    on. ``target`` is a spectrum or spectrum file and ``line_mask`` a mask set or mask file
    naming the element; its windows are measured with ``profile`` (``snr``, ``cont_pct``
    errors) and each is identified as the element's neutral or singly ionised stage from
    the configured line list. ``measured`` continues a :class:`~tumult.measurement.BalanceResult`
    instead, keeping its identity and annotations: its lines, as curated, are solved on
    without measuring again, and a result without lines has ``line_mask`` measured into it
    with its own profile and errors.

    ``method="cog"`` determines per-line abundances by the adaptive curve of growth of
    :func:`tumult.abundance.fit` at every trial star, with the same screens (at the
    starting star) and ``background_lines`` choice. :func:`solve_coupled` steps all free
    parameters through the full diagnostic Jacobian. Every trial star also interpolates its
    rungs at one probe star per free parameter and re-inverts each line there; the
    differences of the diagnostics are the Jacobian, so no probe needs a curve-of-growth
    solve of its own. A probe moves by the largest of :data:`FD_STEPS` (``steps``
    overrides it) and its first :data:`PROBE_HALVINGS` halvings that keeps it on the
    star's own atmosphere nodes, so the Jacobian is the slope of the star's interpolation
    cell. Where the star sits on a grid value, the probe adds the neighbouring nodes: on the
    side its own diagnostic moves it towards (:data:`DIAG_SIGN`) when the two sides need
    different nodes. ``free`` maps parameters to whether they are fitted
    (all by default) and ``tols`` overrides the convergence tolerances.

    Returns a :class:`~tumult.measurement.BalanceResult` with the diagnostics, uncertainties,
    Jacobian, history, ``converged`` and ``termination``; ``request.update(star=result)``
    continues from the solved parameters.
    """
    from copy import copy, deepcopy
    from . import abundance as ab
    from .departures import DepartureGrid, Departures
    from .measurement import (BalanceResult, STAGE_MEASURED, STAGE_SOLVED, lines_key,
                              mask_key)
    from .synthesis.api import NodeStore
    from .synthesis.inputs import AtmosphereGrid
    from .synthesis.operations import quadrature_atmosphere
    from .synthesis.request import SynthRequest
    from .synthesis.workflows import abund_table, node_selector

    if method != "cog":
        raise ValueError(f"Unknown balance abundance method {method!r}; use 'cog'.")
    if not isinstance(request, SynthRequest):
        raise TypeError("sparams.ionbalance expects a SynthRequest")
    free = {k: bool((free or {}).get(k, True)) for k in PARAM_KEYS}
    if any(free.values()) and not isinstance(request.atmosphere, AtmosphereGrid):
        raise ValueError("Stellar parameters can only be fitted on an atmosphere grid; hold "
                         "them all (free={...: False}) to evaluate the diagnostics.")
    if isinstance(request.departures, (Departures, DepartureGrid)):
        raise ValueError("Frozen departures do not follow the star; use departures.request.")
    request.check()
    spectrum, target_record = ab._target_spectrum(target)
    report = progress or (lambda _result: None)
    config = request.config
    mset = None if measured is not None and measured.lines else ab._line_mask_set(line_mask)
    if measured is not None:
        result = deepcopy(measured)
        if mset is not None:
            result.element = mset.element
    else:
        result = BalanceResult(element=mset.element, profile=profile, snr=snr,
                               cont_pct=cont_pct, maskset_name=mset.name,
                               maskset_key=mask_key(mset))
    element, profile = result.element, result.profile
    result.params0, result.free, result.tols = dict(request.star), free, dict(tols or {})
    result.method, result.solver = method, "coupled"
    result.stage, result.history, result.diags, result.sigmas = STAGE_MEASURED, [], {}, {}
    result.converged = False
    result.request, result.target = request, target_record
    result.cog_settings = {"background_lines": bool(background_lines), "limit": float(limit)}
    if mset is not None:
        measure_balance_lines(result, [m for m in mset.masks if m.include], spectrum,
                              str(config.resolve(config.linelist)) if config.linelist else None,
                              atmosphere=quadrature_atmosphere(request))
    report(result)

    solar = abund_table(config.to_prefs())
    solar_x = solar.get(element)
    if solar_x is None:
        raise ValueError(f"No solar {element} abundance; configure config.abundances.")
    result.solar_x = float(solar_x)
    result.termination, result.converged = "running", False

    # one node store for the whole solve: a trial star reuses every rung run an earlier
    # one computed on the atmosphere nodes they share
    store = NodeStore()

    bounds = grid_bounds(request.atmosphere.record(config))
    step_sizes = {**FD_STEPS, **(steps or {})}

    def context_at(star, lines, probes=()):
        return ab._LatticeContext(request.update(star=star), lines, element,
                                  background_lines=background_lines, limit=limit,
                                  cancel_token=cancel_token, store=store, probes=probes)

    ab.screen_lines(context_at(request.star, result.lines), blend_max=blend_max,
                    weak_frac=weak_frac)
    screened = deepcopy(result.lines)
    report(result)

    def at_probe(lines, probe, found):
        """The diagnostics with every solved line at its probe abundance and error (the
        diagnostics' weights move with the star too), or None when a line leaves its curve
        of growth there."""
        moved = []
        for ln in lines:
            a, sigma = found.get(id(ln), (np.nan, np.nan))
            ln = copy(ln)
            if ln.include and np.isfinite(ln.a):
                if not np.isfinite(a):
                    return None
                ln.a, ln.sigma_a = a, sigma
            moved.append(ln)
        return diagnostics(moved, probe["feh"], solar_x)

    select = node_selector(request, store)

    def shifted(params, key, step):
        value = params[key] + step
        lo, hi = bounds.get(key, (-np.inf, np.inf))
        return {**params, key: float(value)} if lo <= value <= hi else None

    def within(params, key, sign, allowed):
        """The farthest probe along ``sign`` (at most one step) whose nodes are all in
        ``allowed``, or None."""
        for n in range(PROBE_HALVINGS + 1):
            probe = shifted(params, key, sign * step_sizes[key] * 0.5 ** n)
            if probe is not None and select(probe) <= allowed:
                return probe
        return None

    def plan_probes(params):
        """Probes on the star's own nodes where it lies inside a cell, and for the axes it
        is pinned on whose two sides need different nodes, the nodes each side reaches:
        ``(probes, {key: {sign: nodes}})``."""
        own = select(params)
        probes, pinned = {}, {}
        for k in PARAM_KEYS:
            if not free[k]:
                continue
            probe = within(params, k, 1, own) or within(params, k, -1, own)
            if probe is not None:
                probes[k] = probe
                continue
            tiny = step_sizes[k] * 0.5 ** PROBE_HALVINGS
            reach = {sign: own | select(probe) for sign in (1, -1)
                     if (probe := shifted(params, k, sign * tiny)) is not None}
            if len(reach) == 2 and reach[1] != reach[-1]:
                pinned[k] = reach
            elif reach:
                # one side left, or both sides reach the same nodes
                sign = 1 if 1 in reach else -1
                probes[k] = within(params, k, sign, reach[sign])
        return probes, pinned

    def solve(params, probes):
        """The star's solved lines, diagnostics and Jacobian columns from ``probes``."""
        lines = deepcopy(screened)
        context = context_at(params, lines, probes=list(probes.values()))
        ab.fit_lines_cog(context, limit=limit, max_rounds=max_rounds, profile=profile)
        # [X/Fe] -> A(X); a probe converts at its own metallicity
        found = [{i: (xfe_to_abund(x, element, probe["feh"], solar), sigma)
                  for i, (x, sigma) in solutions.items()}
                 for probe, solutions in zip(probes.values(), context.probe_solutions)]
        for ln in lines:
            if ln.include and np.isfinite(ln.a):
                shift = params["feh"] + solar_x
                ln.a = xfe_to_abund(ln.a, element, params["feh"], solar)
                ln.a_grid = tuple(a + shift for a in ln.a_grid)
                ln.cog_specs = {a + shift: spec for a, spec in ln.cog_specs.items()}
                ln.cog_origin += shift
        diags = diagnostics(lines, params["feh"], solar_x)
        columns = {}
        for (k, probe), solutions in zip(probes.items(), found):
            moved = at_probe(lines, probe, solutions)
            if moved is not None:
                columns[k] = {d: (moved[d] - diags[d]) / (probe[k] - params[k])
                              for d in DIAG_KEYS}
        return lines, diags, columns

    def evaluate(params):
        probes, pinned = plan_probes(params)
        lines, diags, columns = solve(params, probes)
        if pinned:
            # a star on a grid value along an axis whose two sides need different nodes:
            # probe on the side its own diagnostic says the first step goes
            extra = {}
            for k, reach in pinned.items():
                residual = diags.get(DIAG_FOR[k], np.nan)
                sign = -1 if np.isfinite(residual) and residual * DIAG_SIGN[k] > 0 else 1
                extra[k] = within(params, k, sign, reach[sign])
            columns.update(solve(params, extra)[2])
        jacobian = {d: {k: column[d] for k, column in columns.items()} for d in DIAG_KEYS}
        return diags, lines, jacobian

    def on_evaluation(history):
        result.history = [dict(item) for item in history]
        report(result)

    try:
        outcome = solve_coupled(
            evaluate, dict(request.star), free=free, bounds=bounds, tols=tols,
            max_iter=max_iter, steps=steps, report=on_evaluation,
            cancel_token=cancel_token)
    except Exception as exc:
        from .synthesis.cancel import CancelledError
        result.converged = False
        result.termination = "cancelled" if isinstance(exc, CancelledError) else "failed"
        raise
    result.lines = outcome.state
    result.params, result.diags = outcome.params, outcome.diags
    result.history = outcome.history
    result.jacobian, result.covariance = outcome.jacobian, outcome.covariance
    result.sigmas = outcome.sigmas
    result.termination, result.converged = outcome.termination, outcome.converged
    result.solved_key = lines_key(result.lines)
    result.stage = STAGE_SOLVED
    result.name = result.label()
    report(result)
    return result
