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

The solver is a damped quasi-Newton fixed point on the four decoupled (diagnostic,
parameter) pairs: iteration 1 uses conservative default sensitivities, later iterations
refine each diagonal Jacobian entry by secant from the observed diagnostic change
(sign-guarded and magnitude-clamped, since the diagonal secant is contaminated by
cross-coupling). Steps are clamped per iteration and to the atmosphere grid's hull.

A parameter may be held *fixed*: it then never steps, and its diagnostic is still computed
and reported but does not gate convergence — the residual at a fixed parameter is the
interesting output ("at your log g the ionization balance is off by 0.08 dex"). Fixing all
four turns the solve into a single diagnostic evaluation.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from .abundance import reduced_ew

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
    """True when every finite gating diagnostic is inside its tolerance.

    An all-NaN set is not converged (nothing was measured), and neither is a solve with no
    free parameter — that is a diagnostic evaluation, and the caller labels it as one.
    """
    tols = {**DEFAULT_TOLS, **(tols or {})}
    keys = free_diags(free)
    finite = [k for k in keys if np.isfinite(diags.get(k, float("nan")))]
    if not finite:
        return False
    return all(abs(diags[k]) <= tols[k] for k in finite)


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


def identify_catalog(star, lo, hi, linelist):
    """``(frame, note)``: the transitions in ``[lo, hi]`` Å ranked for ``star``.

    HOP strengths on the reference atmosphere are what makes the identification a statement
    about *this* star; the master line list is the fallback when HOP cannot be run (no line
    list configured, or the band fails), and then the note says so, because a log gf-ordered
    identification is a much weaker claim.
    """
    from . import linestrength
    from .linelist import load_master
    try:
        frame = linestrength.compute_band(star, str(linelist),
                                          *linestrength.band_bounds(lo, hi))
        if frame is not None and len(frame):
            return frame, f"HOP strengths on {star.label}"
    except Exception:
        logging.exception("Line-strength band failed; falling back to the line list")
    try:
        return load_master(linelist), "line-list order (no HOP strengths)"
    except Exception:
        logging.exception("Master line list unreadable")
        return None, "no line list — identities unknown"


def measure_balance_lines(result, masks, spec, linelist):
    """Fill ``result.lines`` from the windows, identified on the entry point's atmosphere.

    Raises :class:`ValueError` when no window could be identified as either stage of the
    element: the balance reads the *stage* off this identification, and the mask set's
    element cannot supply one. Windows that are dominated by another species, or whose
    dominant transition is the third stage or above, are switched off rather than dropped —
    the flag stays visible and the caller's toggle overrides it.
    """
    from . import linediag
    star = entry_star(result.params0)
    lo = min(m.lmin for m in masks)
    hi = max(m.lmax for m in masks)
    catalog, note = identify_catalog(star, lo, hi, linelist)
    linediag.measure_masks(spec, masks, profile=result.profile, snr=result.snr,
                           cont_pct=result.cont_pct, catalog=catalog,
                           element=result.element)
    result.ident_star, result.ident_note = star.label, note

    lines = []
    for m in masks:
        ln = linediag.line_from_mask(m)
        stage = is_ion_stage(ln.ion, 1) or is_ion_stage(ln.ion, 2)
        if "off_element" in ln.flags or not stage:
            ln.include = False
            if not stage and "off_element" not in ln.flags:
                ln.flags = ln.flags + ("no_stage",)
        lines.append(ln)
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
    from .synthesis.workflows import DEFAULT_PATTERN, discover_atmos_grid
    if atmos_file or not grid.get("folder"):
        return {"vmic": VMIC_BOUNDS}
    try:
        files = discover_atmos_grid(grid["folder"], grid.get("pattern") or DEFAULT_PATTERN)
    except Exception:
        return {"vmic": VMIC_BOUNDS}
    keys = np.array([k[:3] for k in files], dtype=float)
    bounds = {name: (float(keys[:, i].min()), float(keys[:, i].max()))
              for i, name in enumerate(("teff", "logg", "feh"))}
    vmics = [k[3] for k in files if len(k) > 3]
    bounds["vmic"] = (float(min(vmics)), float(max(vmics))) if vmics else VMIC_BOUNDS
    return bounds


def _strip_element(text, element):
    """Drop any user ``X=`` term from the base abundance overrides (the solver owns it)."""
    pre = f"{str(element).lower()}="
    return " ".join(t for t in str(text or "").split() if not t.lower().startswith(pre))


def solve_ionisation_balance(request, result, opts=None, *, progress=None, cancel_token=None):
    """Iterate the excitation-ionization balance to convergence on ``result.lines``.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest`; ``result`` a
    :class:`~tumult.measurement.Measurement` already carrying measured lines, an entry point
    (``params0``), which parameters are free and the tolerances. ``opts`` sets the ladder
    (``span``, ``points``) and ``max_iter``.

    Each iteration synthesises a per-line abundance ladder at the trial parameters, inverts
    the measured EW through it, and steps the free parameters on the resulting diagnostics.
    ``progress(result)`` is called as each line and each iteration lands, so a caller can
    stream partial results; it must not raise.
    """
    from . import abundance as ab
    from .abundance import SYNTH_PAD_AA, xfe_to_abund
    from .synthesis import workflows as synth
    from .synthesis.cancel import CancelledError

    opts = {**DEFAULT_SOLVE_OPTS, **(opts or {})}
    report = progress or (lambda _result: None)
    prefs, grid, atmos_file = request.prefs, request.grid, request.atmos_file
    element = result.element
    base_rv = dict(request.run_vals)
    relative = bool(base_rv.get("relative_abundances", False))
    base_abund = _strip_element(base_rv.get("abundances", ""), element)
    solar = synth.abund_table(prefs)
    solar_x = solar.get(element)
    if solar_x is None:
        raise ValueError(f"No solar {element} abundance — configure the abundance file.")

    free = dict(result.free) or {k: True for k in PARAM_KEYS}
    any_free = any(free.values())
    params = dict(result.params0)
    state = SolverState(params=dict(params), bounds=grid_bounds(grid, atmos_file))
    result.solar_x = float(solar_x)

    def rv_point(a_rel, feh, only=True, exclude=False):
        rv = dict(base_rv)
        rv.pop("only_elements", None)
        rv.pop("exclude_elements", None)
        if exclude:
            rv["abundances"] = base_abund
            rv["exclude_elements"] = element
        else:
            val = a_rel if relative else xfe_to_abund(a_rel, element, feh, solar)
            rv["abundances"] = " ".join(
                t for t in (base_abund, f"{element}={val:.4f}") if t)
            if only:
                rv["only_elements"] = element
        rv["relative_abundances"] = relative
        return rv

    def batch(rv_list, line, target):
        return synth.compute_batch(
            target, prefs, rv_list, line.window[0] - SYNTH_PAD_AA,
            line.window[1] + SYNTH_PAD_AA, grid=grid, atmos_file=atmos_file,
            cancel_token=cancel_token)

    # each line's ladder centre carries over between iterations ([Fe/H]-relative)
    centres = {}
    diags = {}
    # with every parameter held there is nothing to step: one pass produces the diagnostics
    # at the given parameters, which is the point of holding them.
    max_iter = int(opts["max_iter"]) if any_free else 1
    span, npts = float(opts["span"]), int(opts["points"])
    for it in range(1, max_iter + 1):
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        target = tuple(params[k] for k in PARAM_KEYS)
        feh = params["feh"]
        for i, ln in enumerate(result.lines):
            if cancel_token is not None and cancel_token.cancelled:
                raise CancelledError()
            if not ln.include or not np.isfinite(ln.ew) or ln.ew <= 0:
                continue
            c = centres.get(i, 0.0)
            ladder = np.linspace(c - span, c + span, npts)
            screen = it == 1        # blend screen rides along on iteration 1
            try:
                # the ladder goes in on its own: its points differ in nothing but the
                # element's abundance, which is what lets compute_spectra fold them into a
                # single run's set axis. The blend screen's two syntheses switch elements on
                # and off, so mixing them into this batch would unfold the ladder.
                specs = batch([rv_point(a, feh) for a in ladder], ln, target)
                if screen:
                    full, excl = batch([rv_point(c, feh, only=False),
                                        rv_point(c, feh, exclude=True)], ln, target)
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
                a_grid, ew_grid = ab.build_cog(
                    None, ladder, window=ln.fit_window, profile=result.profile,
                    x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs: _s)
                ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
                if inv["out_of_range"]:
                    # recentre a doubled ladder on the extrapolated abundance and try once
                    # more — the guess may simply be far off at first
                    c2 = float(np.clip(inv["a"], c - 3 * span, c + 3 * span))
                    ladder2 = np.linspace(c2 - 2 * span, c2 + 2 * span, npts)
                    specs2 = batch([rv_point(a, feh) for a in ladder2], ln, target)
                    a_grid, ew_grid = ab.build_cog(
                        None, ladder2, window=ln.fit_window, profile=result.profile,
                        x0=ln.lam0, batch_synth_fn=lambda _a, _s=specs2: _s)
                    inv = ab.cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
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
        if not any(np.isfinite(diags.get(k, float("nan"))) for k in DIAG_KEYS):
            # not one line could be inverted, so every diagnostic is NaN and every parameter
            # holds — iterating again would synthesise the same thing and learn the same
            # nothing. Usually an entry point too far from the star for the curve-of-growth
            # ladder to reach the measured EWs.
            break
        if any_free and converged(diags, result.tols, free):
            result.converged = True
            break
        if not any_free:
            break
        params = step(state, diags, free=free)

    result.jac = dict(state.jac)
    result.sigmas = uncertainties(diags, state.jac, free)
    return result
