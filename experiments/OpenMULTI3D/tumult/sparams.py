"""Excitation–ionization balance numerics (pure).

The classic spectroscopic determination of (Teff, logg, [Fe/H], vmic) from per-line
abundances of one element in two ionization stages: Teff zeroes the slope of A(X I) vs
excitation energy, vmic the slope of A(X I) vs reduced EW, logg the X I/X II ionization
offset, and the metallicity the consistency of the mean A(X) with the model's. Everything
in the first half is numpy-only; :func:`solve_ionisation_balance` at the end drives the
syntheses that turn measured EWs into the per-line abundances those numerics consume.

Nothing here is Fe-specific. Fe is the element the method is *usually* run on, and its
sensitivities are the ones :data:`DEFAULT_JAC` seeds every element with (they only have to
point the first step the right way; the secant refinement takes over from iteration 2).
For an element other than Fe the metallicity axis closes only under a scaled-solar
assumption — ``d_feh`` compares ⟨A(X)⟩ with ``A_solar(X) + [Fe/H]_model``, so the caller is
asserting [X/Fe] = 0. That is exactly right for Fe, defensible for most iron-peak elements
and known to be wrong for the α elements at low metallicity; the choice belongs to the
caller, which is why the solar abundance arrives as an argument.

Two solvers step the parameters. The *diagonal* solver (:func:`solve_ionisation_balance`)
is a damped quasi-Newton fixed point on the four decoupled (diagnostic, parameter) pairs:
iteration 1 uses conservative default sensitivities, later iterations refine each diagonal
Jacobian entry by secant from the observed diagnostic change (sign-guarded and
magnitude-clamped, since the diagonal secant is contaminated by cross-coupling). The
*coupled* solver (:func:`solve_coupled`, used by :func:`ionbalance`) treats the free
diagnostics as one residual vector with a full Jacobian, so a Teff step that moves the
ionisation offset is accounted for. Both clamp steps per iteration and to the atmosphere
grid's hull.

A parameter may be held *fixed*: it then never steps, and its diagnostic is still computed
and reported but does not gate convergence — the residual at a fixed parameter is the
interesting output ("at your log g the ionization balance is off by 0.08 dex"). Fixing all
four turns the solve into a single diagnostic evaluation.
"""
from __future__ import annotations

import logging
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from .abundance import reduced_ew, xfe_to_abund

PARAM_KEYS = ("teff", "logg", "feh", "vmic")
DIAG_KEYS = ("s_exc", "s_rew", "d_ion", "d_feh")

# which diagnostic drives which parameter (the classic decoupling)
DIAG_FOR = {"teff": "s_exc", "vmic": "s_rew", "logg": "d_ion", "feh": "d_feh"}

# conservative default sensitivities d(diagnostic)/d(parameter) for an FGK star, used
# whatever the element: they only have to point iteration 1 the right way, and the secant
# refinement below owns every step after that (widen JAC_CLAMP if some element's first step
# ever walks the wrong way). The
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


def step(state, diags, clamps=None, free=None):
    """One damped quasi-Newton update; returns the new params dict (and mutates ``state``).

    Each free parameter moves by ``-diag / J`` along its own diagnostic, clamped per
    iteration and to the grid bounds. Before stepping, each Jacobian entry is refined by
    secant from the previous iteration's (param, diagnostic) pair when that parameter
    actually moved — guarded to the default's sign and a magnitude window, because the
    diagonal secant absorbs cross-coupling between the axes. A NaN diagnostic (e.g. no
    second-stage lines survived) holds its parameter, and so does ``free[key] is False``.
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
        if free is not None and not free.get(key, True):
            continue
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


def uncertainties(diags, jac, free=None):
    """Parameter 1σ from the converged diagnostics: ``σ_p = σ_diag / |J_pp|``.

    First-order and diagonal-only — the honest formal error of the balance method
    itself, not a full covariance. NaN where the diagnostic error (or Jacobian) is
    unusable, and NaN for a parameter that was held fixed: it was not determined here, so
    it has no error of its own (what its own diagnostic says about it is the residual).
    """
    out = {}
    for key in PARAM_KEYS:
        if free is not None and not free.get(key, True):
            out[key] = float("nan")
            continue
        dk = DIAG_FOR[key]
        err = diags.get(f"{dk}_err", float("nan"))
        j = jac.get(key, float("nan"))
        out[key] = float(abs(err / j)) if np.isfinite(err) and np.isfinite(j) and j != 0 \
            else float("nan")
    return out


# --- the synthesis-driven solve ------------------------------------------------------------
# Everything above is numpy on measured lines. What follows drives a synthesis: it inverts
# each line's EW to an abundance through a curve of growth computed at the trial parameters,
# turns those abundances into the four diagnostics above, and steps. It is here rather than
# in the GUI because nothing in it is about rendering — a script, a notebook or an MCP
# adapter runs the same solve by handing it a request and a measurement.

BLEND_MAX = 0.30           # iteration-1 blend screen: contaminant EW fraction ceiling
VMIC_BOUNDS = (0.1, 5.0)   # vmic is a synthesis input, not a grid axis — keep it sane

DEFAULT_SOLVE_OPTS = {"max_iter": 8, "span": 0.4, "points": 3, "tols": None}


def _is_shutdown(exc):
    """True for a pool refusing new work because the interpreter is tearing down.

    A background solve outliving its caller sees ``cannot schedule new futures after
    interpreter shutdown``. That is a cancel, not a per-line synthesis failure, so the loop
    aborts instead of re-issuing (and re-failing) every remaining line.
    """
    return isinstance(exc, RuntimeError) and "interpreter shutdown" in str(exc)


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


def solve_ionisation_balance(request, result, opts=None, *, run_values=None, progress=None,
                             cancel_token=None):
    """Iterate the excitation-ionization balance to convergence on ``result.lines``.

    ``request`` is the :class:`~tumult.synthesis.SynthRequest` whose star each iteration
    replaces and ``run_values`` what its syntheses add (:func:`tumult.synthesis.run_values`);
    ``result`` a :class:`~tumult.measurement.BalanceResult` already carrying measured lines,
    an entry point (``params0``), which parameters are free and the tolerances. ``opts``
    sets the ladder (``span``, ``points``) and ``max_iter``.

    Each iteration synthesises a per-line abundance ladder at the trial parameters, inverts
    the measured EW through it, and steps the free parameters on the resulting diagnostics.
    ``progress(result)`` is called as each line and each iteration lands, so a caller can
    stream partial results; it must not raise.
    """
    from .synthesis.cancel import CancelledError
    try:
        solved = _solve_balance(request, result, opts, run_values, progress, cancel_token)
    except CancelledError:
        result.termination, result.converged = "cancelled", False
        raise
    except Exception:
        result.termination, result.converged = "failed", False
        raise
    result.stage = "solved"
    return solved


def _solve_balance(request, result, opts, run_values, progress, cancel_token):
    from . import abundance as ab
    from .abundance import abundance_override, padded_window, strip_abundance, xfe_to_abund
    from .synthesis import workflows as synth
    from .synthesis.cancel import CancelledError
    from .synthesis.inputs import AtmosphereGrid
    from .synthesis.operations import run_values as complete_values

    opts = {**DEFAULT_SOLVE_OPTS, **(opts or {})}
    report = progress or (lambda _result: None)
    config = request.config
    prefs = config.to_prefs()
    grid = request.atmosphere.record(config) \
        if isinstance(request.atmosphere, AtmosphereGrid) else {}
    element = result.element
    extra = dict(run_values or {})
    base_rv = complete_values(request, **extra)
    relative = bool(base_rv.get("relative_abundances", True))
    base_abund = strip_abundance(base_rv.get("abundances", ""), element)
    solar = synth.abund_table(prefs)
    solar_x = solar.get(element)
    if solar_x is None:
        raise ValueError(f"No solar {element} abundance — configure the abundance file.")

    free = dict(result.free) or {k: True for k in PARAM_KEYS}
    any_free = any(free.values())
    params = dict(result.params0)
    state = SolverState(params=dict(params), bounds=grid_bounds(grid))
    result.request = request.update(star=params)
    result.solar_x = float(solar_x)
    result.converged = False
    result.termination = "running"
    result.history = []
    result.params = dict(params)
    result.diags = {}
    if int(opts["max_iter"]) < 1 or float(opts["span"]) <= 0 or int(opts["points"]) < 2:
        raise ValueError("Balance needs max_iter >= 1, span > 0 and points >= 2.")

    # a request with ``only_elements`` (background lines switched off) synthesises the
    # curve-of-growth ladder from the element's lines alone; the blend screen always compares
    # the full synthesis with the element excluded
    ladder_only = bool(base_rv.get("only_elements"))

    def rv_point(a_rel, feh, only=None, exclude=False, save_linelist=True):
        only = ladder_only if only is None else only
        rv = dict(extra)
        # linelist_out.txt lists the transitions the run kept, with strengths at its own
        # composition: a diagnostic, never an input. Membership is the same across a
        # ladder for every line the balance uses, so one rung per node keeps it and the
        # rest pass False rather than rewrite ~150 kB each. Everything else keeps it.
        rv["save_linelist"] = bool(save_linelist)
        rv["only_elements"] = None
        rv["exclude_elements"] = None
        if exclude:
            rv["abundances"] = base_abund
            rv["exclude_elements"] = element
        else:
            rv["abundances"] = abundance_override(
                element, a_rel, base=base_abund, relative=relative, feh=feh, solar=solar)
            if only:
                rv["only_elements"] = element
        rv["relative_abundances"] = relative
        return rv

    def batch(rv_list, line, point):
        return synth.compute_batch(point, rv_list, *padded_window(line.window),
                                   cancel_token=cancel_token)

    def ladder_specs(ladder, centre, feh, lines, point, *, crop=None):
        """One synthesis per lattice rung, each covering every line's window.

        A rung is its own run so that its identity is (atmosphere node, abundance,
        windows) and nothing else: the same rung recurs at the next iteration and is
        read from the cache instead of resynthesised. Folding the ladder onto one run's
        set axis instead would amortise the continuum and LTE sets across the rungs,
        but its cache key is the entire ladder, which no later iteration reproduces.

        The rungs still go out together (``fold=False``): every (rung x node) run is
        issued in one batch, so the pool sees the whole iteration's work at once rather
        than one rung's nodes at a time.

        ``crop`` describes the ladder's strongest rung, so every rung on a node reads the
        same cropped line list and no rung waits for another. Only the rung nearest
        scaled-solar keeps its transition dump, which is the same diagnostic on all of them.
        """
        windows = [padded_window(ln.window) for ln in lines]
        rungs = [rv_point(a, feh, save_linelist=(a == centre)) for a in ladder]
        return synth.compute_window_batch(point, rungs, windows, cancel_token=cancel_token,
                                          crop=crop, fold=False)

    def shared_ladder(active, centres_now):
        """The lattice rungs covering every active line's bracket.

        Rungs sit on a global lattice of ``step`` dex, so a rung is the same run
        whenever it is asked for again and the ladder need not repeat between
        iterations for its points to be reused. The range spans the line-to-line
        scatter plus each line's own +/-span, with both edges snapped outward; that
        scatter is the balance signal, so it is coverage rather than waste. No
        symmetry is imposed -- the lattice is absolute, not centred on the mean.
        """
        cs = [centres_now.get(i, 0.0) for i, _ln in active]
        lo = np.floor((min(cs, default=0.0) - span) / lattice) * lattice
        hi = np.ceil((max(cs, default=0.0) + span) / lattice) * lattice
        ladder = lo + lattice * np.arange(int(round((hi - lo) / lattice)) + 1)
        # the blend screen rides on the rung nearest scaled-solar, which is also the
        # one rung that writes the transition dump
        centre = float(ladder[np.argmin(np.abs(ladder))])
        return centre, ladder

    # each line's ladder centre carries over between iterations ([Fe/H]-relative)
    centres = {}
    diags = {}
    # with every parameter held there is nothing to step: one pass produces the diagnostics
    # at the given parameters, which is the point of holding them.
    max_iter = int(opts["max_iter"]) if any_free else 1
    span, npts = float(opts["span"]), int(opts["points"])
    # not "step": that name belongs to the parameter-stepping function above
    lattice = float(opts.get("abund_step") or config.abundance_step)
    if lattice <= 0:
        raise ValueError("Balance needs a positive abundance lattice step.")
    for it in range(1, max_iter + 1):
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        point = request.update(star=params)
        feh = params["feh"]
        for ln in result.lines:
            if ln.include:
                ln.a = float("nan")
        active = [(i, ln) for i, ln in enumerate(result.lines)
                  if ln.include and np.isfinite(ln.ew) and ln.ew > 0]
        if not active:
            break
        # One synthesis for the whole iteration: the rungs ride the run's set axis and the
        # lines its output patches. The ladder is shared, so it is centred on the lines'
        # mean and widened to reach the outliers (see shared_ladder).
        c, ladder = shared_ladder(active, centres)
        # every rung on a node reads the same line list, cropped at the top rung: the
        # selection covers the weaker rungs, so it is a property of the node
        crop = {"abundances": abundance_override(element, float(ladder[-1]),
                                                 base=base_abund, relative=relative,
                                                 feh=feh, solar=solar),
                "relative_abundances": relative}
        iter_specs = None
        try:
            iter_specs = ladder_specs(ladder, c, feh, [ln for _i, ln in active], point,
                                      crop=crop)
        except CancelledError:
            raise
        except Exception as exc:
            if (cancel_token is not None and cancel_token.cancelled) or _is_shutdown(exc):
                raise
            # fall back to the per-line path rather than losing the iteration
            logging.exception("shared-ladder synthesis failed; falling back per line")
        for wi, (i, ln) in enumerate(active):
            if cancel_token is not None and cancel_token.cancelled:
                raise CancelledError()
            screen = it == 1        # blend screen rides along on iteration 1
            try:
                if iter_specs is not None:
                    specs = [iter_specs[vi][wi] for vi in range(len(ladder))]
                else:
                    # the ladder goes in on its own: its points differ in nothing but the
                    # element's abundance, which is what lets compute_spectra fold them
                    # into a single run's set axis. The blend screen's two syntheses switch
                    # elements on and off, so mixing them in would unfold the ladder.
                    specs = batch([rv_point(a, feh) for a in ladder], ln, point)
                if screen:
                    full, excl = batch([rv_point(c, feh, only=False),
                                        rv_point(c, feh, exclude=True)], ln, point)
                    try:
                        blend, _, _ = ab.blend_fraction(full, excl, window=ln.fit_window,
                                                        profile=result.profile, x0=ln.lam0)
                    except Exception:
                        # the screen is an extra, not the measurement: a line whose blend
                        # cannot be measured is reported unscreened, not lost
                        logging.exception("blend screen failed at %.3f", ln.lam0)
                        blend = float("nan")
                    ln.blend_frac = float(blend)
                    if np.isfinite(blend) and blend > BLEND_MAX:
                        if "blended" not in ln.flags:
                            ln.flags = ln.flags + ("blended",)
                        ln.include = False
                        report(result)
                        continue
                renorm = bool(getattr(ln.mask, "renorm", False))
                a_grid, ew_grid = ab.build_cog(
                    None, ladder, window=ln.fit_window, profile=result.profile,
                    x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs: _s, renorm=renorm)
                ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
                if inv["out_of_range"]:
                    # recentre a doubled ladder on the extrapolated abundance and try once
                    # more — the guess may simply be far off at first
                    c2 = float(np.clip(inv["a"], c - 3 * span, c + 3 * span))
                    ladder2 = np.linspace(c2 - 2 * span, c2 + 2 * span, npts)
                    # the retry replaces the ladder, so its spectra replace the rungs the
                    # solution spectrum is later read off
                    specs = batch([rv_point(a, feh) for a in ladder2], ln, point)
                    a_grid, ew_grid = ab.build_cog(
                        None, ladder2, window=ln.fit_window, profile=result.profile,
                        x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs: _s, renorm=renorm)
                    inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
                if inv["out_of_range"]:
                    ln.a = float("nan")
                    ln.flags = tuple(f for f in ln.flags if f != "out_of_range") \
                        + ("out_of_range",)
                else:
                    a_rel = float(inv["a"])
                    centres[i] = a_rel
                    ln.a = xfe_to_abund(a_rel, element, feh, solar)
                    ln.sigma_a = float(inv["sigma_a"])
                    ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
                    ln.a_grid, ln.ew_grid = tuple(a_grid), tuple(ew_grid)
                    ln.method = "cog"
                    ln.flags = tuple(f for f in ln.flags if f != "out_of_range")
                    # the line as this iteration's parameters and abundance leave it, read
                    # off the ladder it was inverted through
                    ln.fit_spec = ab.cog_spectrum(a_grid, specs, a_rel)
            except CancelledError:
                raise
            except Exception as exc:
                # a cancel can surface as a non-CancelledError — a killed pool
                # (BrokenProcessPool/RuntimeError) mid-submit, or the executor refusing work
                # once the interpreter is exiting. Abort the whole solve then, rather than
                # logging (and re-issuing) every remaining line.
                if (cancel_token is not None and cancel_token.cancelled) or _is_shutdown(exc):
                    raise CancelledError()
                logging.exception("%s balance: line %.3f failed", element, ln.lam0)
                ln.a = float("nan")
                ln.flags = ("synth_failed",)
            report(result)

        diags = diagnostics(result.lines, feh, solar_x)
        result.history.append({"params": dict(params), "diags": dict(diags)})
        result.params = dict(params)
        result.diags = dict(diags)
        report(result)
        result.termination = balance_status(diags, result.tols, free)
        result.converged = result.termination == "converged"
        if result.termination != "iteration_limit":
            break
        if it < max_iter:
            previous = params
            params = step(state, diags, free=free)
            if params == previous:
                result.termination = "boundary"
                break

    if not result.history:
        result.termination = "undetermined"
    if result.termination == "undetermined" and any(
            ln.include and "synth_failed" in ln.flags for ln in result.lines):
        result.termination = "failed"

    result.jac = dict(state.jac)
    result.sigmas = uncertainties(diags, state.jac, free)
    return result


# --- the coupled solve -----------------------------------------------------------------------
# Numerical contract of solve_coupled:
#
# * residuals: the free diagnostics, each divided by its tolerance, so every component is
#   dimensionless and "converged" means every |r_k| <= 1;
# * variables: each free parameter divided by its finite-difference step (FD_STEPS), so the
#   damping treats one step of every parameter alike;
# * Jacobian: forward differences at the starting point, one evaluation per free parameter,
#   stepping backwards where the forward point leaves the grid; after an accepted step it
#   is updated by Broyden's rank-one formula, and it is recomputed by finite differences
#   where a step is rejected on a Broyden estimate;
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

    ``evaluate(params) -> (diags, state)`` computes the diagnostics at a parameter point;
    ``state`` is whatever the caller wants back for the accepted point (the per-line
    abundances). ``free`` maps parameters to whether they are fitted, ``bounds`` gives
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
        diags, state = evaluate(dict(params))
        history.append({"params": dict(params), "diags": dict(diags), "kind": kind})
        report(history)
        return diags, state

    def residual(diags):
        return np.array([float(diags.get(d, float("nan"))) / tol[d] for d in dkeys])

    def fd_jacobian(params, r):
        """Columns d r / d x at ``params``, x the parameters in units of their step."""
        J = np.empty((len(dkeys), len(keys)))
        for j, k in enumerate(keys):
            probe = dict(params)
            probe[k] = params[k] + scale[j]
            lo, hi = bounds.get(k, (-np.inf, np.inf))
            if probe[k] > hi:
                probe[k] = params[k] - scale[j]
            probe = project(probe)
            dx = (probe[k] - params[k]) / scale[j]
            if dx == 0:
                J[:, j] = np.nan
                continue
            diags, _state = run(probe, "jacobian")
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
    diags, state = run(params, "step")
    if not keys:
        return outcome("diagnostic", None)
    J, fresh, lam = None, False, LM_LAMBDA
    for _it in range(int(max_iter)):
        r = residual(diags)
        if not np.all(np.isfinite(r)):
            return outcome("undetermined", J)
        if converged(diags, tol, free):
            return outcome("converged", J)
        if J is None:
            J, fresh = fd_jacobian(params, r), True
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
        new_diags, new_state = run(trial, "step")
        r_new = residual(new_diags)
        if np.all(np.isfinite(r_new)) and np.linalg.norm(r_new) < np.linalg.norm(r):
            dxm = moved / scale
            J = J + np.outer((r_new - r) - J @ dxm, dxm) / float(dxm @ dxm)
            params, diags, state, fresh = trial, new_diags, new_state, False
            lam = max(lam / 3.0, 1e-8)
        else:
            history[-1]["kind"] = "rejected"
            lam *= 4.0
            if not fresh:
                J = None
    r = residual(diags)
    if not np.all(np.isfinite(r)):
        return outcome("undetermined", J)
    return outcome("converged" if converged(diags, tol, free) else "iteration_limit", J)


def ionbalance(request, *, target, line_mask=None, measured=None, method="cog",
               background_lines=True, solver="coupled", free=None, tols=None,
               max_iter=COUPLED_MAX_ITER, profile="voigt", snr=None, cont_pct=0.0, limit=1.0,
               max_rounds=8, steps=None, blend_max=BLEND_MAX, weak_frac=0.10, progress=None,
               cancel_token=None, **diagonal_opts):
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
    starting star) and ``background_lines`` choice. ``solver="coupled"`` (the default)
    steps all free parameters through the full diagnostic Jacobian (:func:`solve_coupled`);
    ``solver="diagonal"`` is :func:`solve_ionisation_balance` with its own shared ladder
    (``diagonal_opts``: ``span``, ``points``). ``free`` maps parameters to whether they are
    fitted (all by default) and ``tols`` overrides the convergence tolerances.

    Returns a :class:`~tumult.measurement.BalanceResult` with the diagnostics, uncertainties,
    Jacobian, history, ``converged`` and ``termination``; ``request.update(star=result)``
    continues from the solved parameters.
    """
    from copy import deepcopy
    from . import abundance as ab
    from .departures import DepartureGrid, Departures
    from .measurement import BalanceResult, STAGE_MEASURED, STAGE_SOLVED, mask_key
    from .synthesis.inputs import AtmosphereGrid
    from .synthesis.operations import quadrature_atmosphere
    from .synthesis.request import SynthRequest
    from .synthesis.workflows import abund_table

    if method != "cog":
        raise ValueError(f"Unknown balance abundance method {method!r}; use 'cog'.")
    if solver not in ("coupled", "diagonal"):
        raise ValueError("solver must be 'coupled' or 'diagonal'.")
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
    result.method, result.solver = method, solver
    result.stage, result.history, result.diags, result.sigmas = STAGE_MEASURED, [], {}, {}
    result.converged = False
    result.request, result.target = request, target_record
    if mset is not None:
        measure_balance_lines(result, [m for m in mset.masks if m.include], spectrum,
                              str(config.resolve(config.linelist)) if config.linelist else None,
                              atmosphere=quadrature_atmosphere(request))
    report(result)

    run_values = {} if background_lines else {"only_elements": element}
    if solver == "diagonal":
        opts = {"max_iter": max_iter, **diagonal_opts}
        return solve_ionisation_balance(request, result, opts, run_values=run_values,
                                        progress=progress, cancel_token=cancel_token)

    solar = abund_table(config.to_prefs())
    solar_x = solar.get(element)
    if solar_x is None:
        raise ValueError(f"No solar {element} abundance; configure config.abundances.")
    result.solar_x = float(solar_x)
    result.termination, result.converged = "running", False

    def context_at(star, lines):
        return ab._CogContext(request.update(star=star), lines, element,
                              background_lines=background_lines, limit=limit,
                              cancel_token=cancel_token)

    ab.screen_lines(context_at(request.star, result.lines), blend_max=blend_max,
                    weak_frac=weak_frac)
    screened = deepcopy(result.lines)
    report(result)

    def evaluate(params):
        lines = deepcopy(screened)
        context = context_at(params, lines)
        ab.fit_lines_cog(context, limit=limit, max_rounds=max_rounds, profile=profile)
        for ln in lines:
            if ln.include and np.isfinite(ln.a):
                ln.a = xfe_to_abund(ln.a, element, params["feh"], solar)
        return diagnostics(lines, params["feh"], solar_x), lines

    def on_evaluation(history):
        result.history = [dict(item) for item in history]
        report(result)

    try:
        outcome = solve_coupled(
            evaluate, dict(request.star), free=free, bounds=grid_bounds(request.atmosphere.record(config)),
            tols=tols, max_iter=max_iter, steps=steps, report=on_evaluation,
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
    result.jac = {k: outcome.jacobian[DIAG_FOR[k]][k] for k in PARAM_KEYS
                  if DIAG_FOR[k] in outcome.jacobian}
    result.sigmas = outcome.sigmas
    result.termination, result.converged = outcome.termination, outcome.converged
    result.stage = STAGE_SOLVED
    result.name = result.label()
    report(result)
    return result
