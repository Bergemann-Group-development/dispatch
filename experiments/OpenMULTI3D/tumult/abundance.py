"""Headless, uncertainty-aware abundance and curve-of-growth workflows.

The synthesis-agnostic numerics (``cog_invert``, ``aggregate``, ``reduced_ew``, ``profile_fwhm``)
now live in :mod:`tumult.abundance` and are re-exported here for back-compat. This module keeps the
synthesis-driven, application-facing pieces:

* :func:`measure_line` — stage 1: fit a line on the reference to get its EW, FWHM and windows.
* :func:`fit_line` — stage 2: screen a line (weak / blended), then determine its abundance by
  **directly fitting** the synthetic to the observation (abundance + vshift + vdop) and measure the
  numerical ``dREW/dA`` gradient. Synthesis is injected as callables so this stays engine-agnostic.
* :func:`finalize_line` / :func:`aggregate` — stage 3: fold SNR + continuum-% errors into each
  line's EW, convert to an abundance error through the REW gradient, and combine into an
  inverse-variance element mean.

Error model: SNR noise is independent per pixel and the continuum offset is independent per line, so
both fold into each line's EW error inside ``tumult.utils.get_weq`` (``sigma`` → photon term,
``syst`` → continuum term). The abundance error works in *reduced* EW (``log10 W/λ``, unit
independent): ``σ_REW = σ_W/(W ln10)`` and ``σ_A = σ_REW / |dREW/dA|``.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path

import astropy.units as u
import numpy as np

from .linemask import LineMask
from .spectrum import Spectrum
from ._abundance_numerics import (  # noqa: F401 (public re-exports)
    aggregate,
    aggregate_unweighted,
    abund_to_xfe,
    cog_invert,
    convert_abundances,
    profile_fwhm,
    reduced_ew,
    sigma_a_from_rew,
    velocity_hwhm,
    read_solar_abundances,
    xfe_to_abund,
)
from .measured import MeasuredLines
from .profile import get_profile_class
from .utils import auto_crop_core

# the per-line synthesis / jump window spans SYNTH_HWHM_MULT half-widths (HWHM) either side
# of the *fitted* line centre, but is always widened to reach at least FIT_MARGIN_HWHM HWHM
# beyond the auto-crop fit window, so the fit region stays fully enclosed. This window sets
# both the jump view range and what gets synthesised, so the preview matches the synthesis.
SYNTH_HWHM_MULT = 3.0
FIT_MARGIN_HWHM = 1.0

# extra Angstrom synthesised each side of a line's window, so its wings are computed
# rather than extrapolated from the window edge.
SYNTH_PAD_AA = 0.5


def abundance_override(element, value, *, base="", relative, feh, solar):
    """The ``abundances`` namelist override setting ``element`` to ``value``.

    ``value`` is in the caller's convention: [X/Fe] when ``relative``, else an absolute A(X),
    which is what :func:`xfe_to_abund` converts to. ``base`` is any other override the caller
    holds fixed and is kept ahead of the element's term.

    Written to four decimals because this string is part of a run's namelist and therefore of
    its cache identity: two abundances that differ below the printed precision would collide
    on one cached run. Display formatting belongs elsewhere.
    """
    val = value if relative else xfe_to_abund(value, element, feh, solar)
    base = strip_abundance(base, element)
    return " ".join(t for t in (base, f"{element}={val:.4f}") if t)


def strip_abundance(text, element):
    """Remove an element override while preserving all other terms and their precision."""
    prefix = f"{element}=".lower()
    return " ".join(term for term in str(text or "").replace(",", " ").split()
                    if not term.lower().startswith(prefix))


def padded_window(window, pad=SYNTH_PAD_AA):
    """``window`` widened by ``pad`` Å on each side — what actually gets synthesised.

    A line is measured over its window but must be *computed* over more than that, or its
    wings are extrapolated from the window edge instead of calculated. Every synthesis-driven
    per-line workflow pads the same way, so the padding is stated once here.
    """
    return (window[0] - pad, window[1] + pad)


# --- per-line record + element aggregation -----------------------------------------------

@dataclass
class LineResult:
    """One diagnostic line's measurement + derived abundance, filled progressively.

    A line advances through the wizard's stages, so most fields start unset (``NaN``) and are
    populated as each stage runs: ``ew``/``fwhm``/``window`` at line-selection; ``blend_frac``,
    ``a``, ``vshift``/``vdop``, ``drewda`` and ``fit_spec`` at synthesis; ``ew_err``/``sigma_a``
    and its photon/continuum split at uncertainties.

    ``sigma_a_phot`` / ``sigma_a_cont`` are the photon (SNR) and continuum-offset shares of
    ``sigma_a`` (``sigma_a**2 = sigma_a_phot**2 + sigma_a_cont**2``); keeping them lets
    :func:`aggregate` report what dominates the element error budget. ``fit_spec`` is the best-fit
    synthetic line, overlaid on the observation when the row is jumped to.
    """
    lam0: float                 # central wavelength [AA]
    ion: str                    # species label, e.g. "Fe I"
    e_low: float                # lower-level excitation energy [eV]
    loggf: float
    ew: float = float("nan")    # measured equivalent width [maa]
    ew_err: float = float("nan")  # combined 1σ EW error [maa]
    ew_species: float = float("nan")  # query-species-only EW in the fit window [maa] (weak screen)
    fwhm: float = float("nan")  # fitted profile FWHM [AA]
    mask: LineMask = field(default_factory=LineMask)  # fit [lmin,lmax] + synth [smin,smax] windows
    a: float = float("nan")     # abundance (A(X) or [X/Fe], caller's convention)
    vshift: float = float("nan")  # best-fit radial-velocity shift [km/s]
    vdop: float = float("nan")    # best-fit Gaussian (Doppler) broadening [km/s]
    drewda: float = float("nan")  # numerical dREW/dA gradient at the best fit (per dex)
    sigma_a: float = float("nan")
    sigma_a_phot: float = 0.0
    sigma_a_cont: float = 0.0
    blend_frac: float = float("nan")  # contaminant EW fraction in the window (0 = clean)
    a_grid: tuple = ()          # (optional) curve of growth: abundances …
    ew_grid: tuple = ()         # … and the EW each produced (mÅ)
    cog_specs: dict = field(default_factory=dict)  # offset -> rung Spectrum over the synth
                                # window, so the fit window can be re-read after the run
    fit_spec: object = None     # best-fit synthetic Spectrum, for the jump overlay
    fit_func: object = None     # fitted profile f(λ_Å)->flux, for the EW-area hatch on jump
    method: str = "fit"         # "fit" (direct) or "cog"
    flags: tuple = ()
    include: bool = True        # user toggle: count this line in the element mean

    # --- window views: the fit/synth windows now live on the embedded LineMask, but the
    # pipeline (synth.fit_line, measure_ew, the balance / correction tools) still
    # reads them as ``ln.fit_window`` / ``ln.window`` tuples, so expose them read-only here.
    @property
    def fit_window(self) -> tuple:
        return self.mask.fit_window

    @property
    def window(self) -> tuple:
        return self.mask.synth_window

    # --- session persistence (see AbundanceResult / the session snapshot) -------------
    _TUPLE_FIELDS = ("a_grid", "ew_grid", "flags")

    def to_dict(self):
        """A JSON-safe view. ``fit_spec`` is embedded as a Spectrum dict; ``fit_func`` (a live
        callable) is dropped — it only feeds the jump-overlay hatch and is not restored."""
        d = {}
        for f in fields(self):
            if f.name == "fit_func":
                continue
            v = getattr(self, f.name)
            if f.name == "fit_spec":
                v = v.to_dict() if v is not None else None
            elif f.name == "cog_specs":
                v = {str(a): sp.to_dict() for a, sp in v.items()}
            elif f.name == "mask":
                v = v.to_dict()
            elif f.name in self._TUPLE_FIELDS:
                v = list(v)
            d[f.name] = v
        return d

    @classmethod
    def from_dict(cls, d):
        names = {f.name for f in fields(cls)}
        d = {k: v for k, v in (d or {}).items() if k in names}
        for k in cls._TUPLE_FIELDS:
            if d.get(k) is not None:
                d[k] = tuple(d[k])
        if d.get("fit_spec") is not None:
            d["fit_spec"] = Spectrum.from_dict(d["fit_spec"])
        d["cog_specs"] = {float(a): Spectrum.from_dict(sp)
                          for a, sp in (d.get("cog_specs") or {}).items()}
        if d.get("mask") is not None:
            d["mask"] = LineMask.from_dict(d["mask"])
        return cls(**d)


# --- EW measurement + synthesis-driven helpers (inject synth callables to stay decoupled) ----

def _crop(spec, window):
    return spec.crop(xmin=window[0], xmax=window[1]) if window else spec


def measure_profile(spec, window=None, *, profile, x0=None, **kw):
    """Measure a bound :class:`~tumult.profile.Profile` over ``window`` (``auto_crop`` off).

    Thin wrapper over :meth:`tumult.Spectrum.measure`, which owns the window crop and the
    plain-``x0`` unit coercion. ``sigma`` / ``syst`` (kwargs) fold in the photon / continuum
    error terms. The returned Profile exposes the width (``.maa``), the reduced EW (``.rew``)
    and the error split off one interface.
    """
    return spec.measure(window or None, profile=profile, x0=x0, auto_crop=False, **kw)


def _renormalised(spec, renorm):
    """``spec`` divided by its own local straight continuum, when ``renorm`` asks for it.

    The fit spans the whole spectrum handed in — for a line that is the synth window, wider
    than the window the width is integrated over — so the shoulders outside the feature are
    what sets the slope. The observation gets the same treatment from its mask's markers
    (:meth:`tumult.linemask.LineMask.renorm_knots`), fitted on its own flux rather than
    imported from here: the two sides share the definition of a local continuum, not a level.
    """
    if not renorm:
        return spec
    from .continuum import divide_local
    return divide_local(spec)


def _ew(spec, window=None, profile="trapz", x0=None, renorm=False):
    """Noiseless *width* EW of a (synthetic) spectrum over ``window``, in mÅ."""
    return measure_profile(_renormalised(spec, renorm), window,
                           profile=profile, x0=x0).maa


def _rew(spec, window=None, profile="trapz", x0=None, renorm=False):
    """Noiseless reduced EW ``log10(W/λ)`` of a (synthetic) spectrum over ``window`` (dimensionless)."""
    return measure_profile(_renormalised(spec, renorm), window,
                           profile=profile, x0=x0).rew


def measure_ew(spec, window=None, snr=None, cont_pct=0.0, profile="trapz", x0=None):
    """Measure a line's EW as a :class:`~tumult.profile.Profile` with a photon / continuum error split.

    ``snr`` None (a synthetic reference) drops the photon term; ``cont_pct`` 0 drops the continuum
    term. A single measurement carries both shares — ``.err_phot_maa`` (SNR) and ``.err_cont_maa``
    (continuum offset), combining in quadrature to ``.width_err_maa`` — so no second pass is needed.
    (The REW conversion for the abundance error happens in :func:`sigma_a_from_rew`.)
    """
    sub = _crop(spec, window)
    n = len(sub)
    sigma = np.full(n, 1.0 / float(snr)) if snr else None
    syst = (float(cont_pct) / 100.0) or None
    try:
        return measure_profile(sub, profile=profile, x0=x0, sigma=sigma, syst=syst)
    except ValueError:
        # the fitted profile can fail on a narrow / noisy window; a trapz width is robust
        return measure_profile(sub, profile="trapz", x0=x0, sigma=sigma, syst=syst)


@dataclass(frozen=True)
class EwRefresh:
    """What :func:`refresh_ews` did to a list of lines, by position in that list.

    ``changed`` holds the indices whose ``ew``/``ew_err`` were reassigned, ``skipped`` those
    without a finite fit window (no EW to refresh) and ``failed`` maps the indices whose
    measurement raised to the error message; those lines keep their previous values.
    """
    changed: tuple = ()
    skipped: tuple = ()
    failed: dict = field(default_factory=dict)


def refresh_ews(spec, lines, *, profile, snr=None, cont_pct=0.0):
    """Re-measure each line's EW on ``spec`` over its stored fit window, in place.

    The remeasurement the abundance and ionisation-balance tools share when the
    observation they were measured on changes appearance (the canvas filter). Only ``ew``
    and ``ew_err`` are written; windows, flags, ``include`` and every fitted quantity are
    left as curated. Whether the lines *should* be refreshed (the result's stage, its
    selection) is the caller's decision. A line whose measurement raises is reported in
    the returned :class:`EwRefresh` rather than silently skipped.
    """
    changed, skipped, failed = [], [], {}
    for i, ln in enumerate(lines):
        if not np.isfinite(ln.fit_window[0]):
            skipped.append(i)
            continue
        try:
            prof = measure_ew(spec, window=ln.fit_window, snr=snr, cont_pct=cont_pct,
                              profile=profile, x0=ln.lam0)
        except Exception as exc:
            logging.warning("EW refresh failed for line %d (%.3f A): %s", i, ln.lam0, exc)
            failed[i] = f"{type(exc).__name__}: {exc}"
            continue
        ln.ew, ln.ew_err = prof.maa, prof.width_err_maa
        changed.append(i)
    return EwRefresh(tuple(changed), tuple(skipped), failed)


def build_cog(synth_line_fn, a_ladder, window=None, profile="trapz", x0=None,
              batch_synth_fn=None, renorm=False):
    """Synthesise the line at each abundance in ``a_ladder`` and measure its *width* EW.

    ``synth_line_fn(a) -> Spectrum`` returns the continuum-normalised synthetic line at
    abundance ``a`` (caller's convention — [X/Fe] or A(X)). Returns ``(a_grid, ew_grid)`` (mÅ).

    ``batch_synth_fn(a_list) -> [Spectrum]`` (optional): when given, the whole ladder is
    synthesised in one call so the front-end can fan every point's run across a single process
    pool; the returned spectra must align with ``a_list``. Equivalent to the serial
    ``synth_line_fn`` loop (same per-point spectra), just scheduled concurrently — used when
    omitted.

    ``renorm`` divides each rung by its own local straight continuum before integrating, to
    match an observation measured against its mask's local continuum. Set it from the mask,
    so the ladder and the measurement it is inverted against agree.
    """
    a = np.asarray(a_ladder, dtype=float)
    if batch_synth_fn is not None:
        specs = batch_synth_fn(list(a))
        ews = np.array([_ew(s, window, profile, x0, renorm) for s in specs])
    else:
        ews = np.array([_ew(synth_line_fn(av), window, profile, x0, renorm) for av in a])
    return a, ews


def cog_spectrum(a_ladder, specs, a):
    """The ladder's synthetic line at abundance ``a``, interpolated between its rungs.

    A curve-of-growth solve never synthesises the abundance it lands on — it inverts the
    measured EW through the ladder — so the spectrum that goes with the solution has to be
    read back off the rungs. The two bracketing rungs are combined linearly in flux at the
    inverted abundance, which is what makes a review overlay possible without a synthesis per
    line per iteration; a request outside the ladder returns its nearest rung.

    Returns ``None`` when there is nothing to interpolate (no rungs, a non-finite ``a``, or
    rungs that do not share a wavelength grid).
    """
    a_grid = np.asarray(a_ladder, dtype=float)
    if a_grid.size == 0 or not specs or not np.isfinite(a):
        return None
    order = np.argsort(a_grid)
    a_grid, specs = a_grid[order], [specs[i] for i in order]
    j = int(np.clip(np.searchsorted(a_grid, float(a)), 1, len(a_grid) - 1))
    lo, hi = specs[j - 1], specs[j]
    x_lo, y_lo = lo.values
    x_hi, y_hi = hi.values
    if len(x_lo) != len(x_hi):
        return lo if abs(a - a_grid[j - 1]) <= abs(a - a_grid[j]) else hi
    span = a_grid[j] - a_grid[j - 1]
    w = float(np.clip((float(a) - a_grid[j - 1]) / span, 0.0, 1.0)) if span else 0.0
    y = np.asarray(y_lo, dtype=float) * (1.0 - w) + np.asarray(y_hi, dtype=float) * w
    return Spectrum(np.asarray(x_lo, dtype=float), y, xunit=lo.xunit, yunit=lo.yunit)


def blend_fraction(spec_all, spec_excl, window=None, profile="trapz", x0=None):
    """Contaminant EW fraction in the window: ``EW(element excluded) / EW(all lines)``.

    ``spec_all`` is the normal synthesis; ``spec_excl`` is the same window synthesised with
    the target element removed (``exclude_elements``). 0 ⇒ clean line, →1 ⇒ fully blended.
    Returns ``(blend_frac, ew_all, ew_excl)`` (widths in mÅ).

    Both arms are **integrated**, whatever ``profile`` the caller measures its lines with
    (the argument is kept so callers can pass their profile through without special-casing).
    A fitted profile is the wrong estimator here: the excluded arm of a clean line is a flat
    window, where a Gaussian fit is unconstrained — it does not fail, it returns an absurdly
    broad, shallow line and an EW of thousands of mÅ, which reads back as a fully blended
    line and silently throws a good one away. Integration is well defined for both a deep
    line and an empty window, and using the same estimator on both arms is what makes the
    ratio mean anything.
    """
    ew_all = _ew(spec_all, window, "trapz", x0)
    ew_excl = _ew(spec_excl, window, "trapz", x0)
    frac = float(ew_excl / ew_all) if ew_all > 0 else float("inf")
    return frac, ew_all, ew_excl


# --- mask-free window + diagnostic screen ------------------------------------------------

def auto_window(spec, lam0, exclude_wing=False, fallback=0.4, resolution=None):
    """Inflection-based fittable window around ``lam0`` (mask replacement).

    Wraps :func:`tumult.utils.auto_crop_core` on the reference spectrum and returns
    ``((xmin, xmax), found)`` — the bounds in the spectrum's wavelength unit, and whether a
    real core was resolved. When the detector cannot resolve a core (e.g. too few/flat
    points) the bounds fall back to ``lam0 ± fallback`` Å and ``found`` is ``False`` so the
    caller can reject the line rather than trust the fallback window.
    """
    xx, yy = spec.values
    try:
        xc, _ = auto_crop_core(xx, yy, lam0, exclude_wing=exclude_wing, resolution=resolution)
        if len(xc) >= 3:
            return (float(np.min(xc)), float(np.max(xc))), True
    except Exception:
        pass
    return (float(lam0 - fallback), float(lam0 + fallback)), False


def local_max_window(spec, lam0, window, tol=1e-4):
    """Narrow ``window`` to the flux maxima nearest ``lam0`` on either side.

    From the sample nearest ``lam0`` the walk first descends to the line's own minimum,
    then climbs outwards on each side and stops where the flux has fallen by ``tol`` below
    the highest value passed; the edge is that highest sample. A side whose flux never
    falls keeps the window's edge. Between two blended lines the stop is the flux maximum
    separating them, so the returned span holds the line's own profile and none of its
    neighbour's core. ``tol`` (in the spectrum's normalised flux) lets the walk cross a
    flat continuum or a saturated flat core without stopping on numerical ripple. Returns
    ``(lo, hi)`` in the spectrum's wavelength unit, never wider than ``window``.
    """
    xx, yy = _crop(spec, window).values
    xx, yy = np.asarray(xx, float), np.asarray(yy, float)
    if xx.size < 3:
        return float(window[0]), float(window[1])
    k0 = int(np.argmin(np.abs(xx - lam0)))

    def climb(step):
        k, last = k0, xx.size - 1 if step > 0 else 0
        while k != last and yy[k + step] < yy[k]:
            k += step                       # down to the core
        top = k
        while k != last and yy[k + step] >= yy[top] - tol:
            k += step
            if yy[k] > yy[top]:
                top = k
        return last if k == last else top

    return float(xx[climb(-1)]), float(xx[climb(+1)])


def overlapping_windows(windows):
    """The keys of ``windows`` (``{key: (lo, hi)}``) whose span overlaps another's."""
    items = sorted(windows.items(), key=lambda kv: kv[1][0])
    hit, reach, prev = set(), -np.inf, None
    for key, (lo, hi) in items:
        if lo < reach:
            hit.add(key)
            hit.add(prev)
        if hi > reach:
            reach, prev = hi, key
    return hit


def screen(ew_maa, blend_frac, lam0_aa, *, ew_min=5.0, logrw_max=-4.8, blend_max=0.30):
    """Flag a candidate line as too weak / saturated / blended; return ``(flags, ok)``.

    ``ok`` is True only when no flag fires — i.e. a clean, unsaturated, sufficiently strong
    diagnostic. Thresholds are exposed so the tool can let the user tune them.
    """
    flags = []
    if ew_maa < ew_min:
        flags.append("weak")
    if reduced_ew(ew_maa, lam0_aa) > logrw_max:
        flags.append("strong")
    if blend_frac > blend_max:
        flags.append("blended")
    return tuple(flags), not flags


# --- per-line + per-element orchestration ------------------------------------------------

@dataclass
class AbundanceResult(MeasuredLines):
    """One completed element run: the per-line measurements plus the aggregated mean.

    A plain data container (no GUI deps) so it can be serialised into a session snapshot and
    re-consumed by a later Fe-balance / stellar-parameter solver. ``feh``/``relative`` record
    the abundance convention so ``a`` values are interpretable.
    """
    symbol: str
    summary: dict
    synthesis: dict = field(default_factory=dict)
    termination: str = "measured"
    stage: str = "lines"
    teff: float = float("nan")
    logg: float = float("nan")
    feh: float = float("nan")
    vmic: float = float("nan")
    relative: bool = True              # True → ``a`` is [X/Fe]; False → A(X)
    uid: str = ""
    name: str = ""
    atom_path: str = ""                # the model atom a correction run was computed from
    quadrature: dict = field(default_factory=dict)  # PreparedLineQuadrature the stages share
    with_linelist: bool = False        # the synthesis line list was included in the runs
    linelist_dump: str = ""            # the reference run's linelist_out.txt, when it was



@dataclass
class Thresholds:
    ew_min: float = 5.0         # mÅ — below this a line is too weak to trust
    logrw_max: float = -4.8     # reduced-EW ceiling — above this the line is saturating
    blend_max: float = 0.30     # contaminant EW fraction ceiling
    weak_frac: float = 0.10     # species-only EW must be ≥ this fraction of the measured EW


# --- staged pipeline (the 3-step wizard drives these one stage at a time) -----------------

def measure_line(ref_spec, line, profile="gauss", fallback=0.4, resolution=None,
                 window=None):
    """Stage 1: fit a line on the reference and return the bound Profile plus its windows.

    Independent of any stellar parameters. With ``window`` given, both the fit and the
    synthesis window are that span (the line's own window from its model atom, say) and no
    auto-crop runs; that is also the only way an integrating profile (``"trapz"``) can be
    used here, since it has no fitted centre or FWHM to derive a window from. Returns
    ``(prof, window, fit_window, flags)``:
    ``prof`` is the fitted :class:`tumult.profile.Profile` — it carries the EW (``.maa``),
    the FWHM, the model ``.func`` (for hatching the integrated area on jump) **and the
    fit-covariance EW error** (``.width_err_maa``; with no per-pixel noise given, the
    covariance is residual-scaled, so this error measures the profile/model mismatch —
    nonzero even on a noiseless synthetic spectrum). ``fit_window`` is the auto-crop
    window the profile was actually fitted over, and ``window`` the wider synthesis /
    jump window: ``SYNTH_HWHM_MULT`` HWHM either side of the *fitted* centre, widened as
    needed to reach ``FIT_MARGIN_HWHM`` HWHM beyond ``fit_window`` on each side.

    ``flags`` carries any sanity-check rejection reasons: ``"no_window"`` when
    :func:`auto_window` could not resolve a core — in that case no EW is measured at all
    and ``(None, (NaN, NaN), (NaN, NaN), ("no_window",))`` is returned — and
    ``"off_center"`` when the fitted line centre falls outside the window. A non-empty
    ``flags`` means the caller should keep the line but mark it rejected. Raises only on
    a hard fit failure (too few points), so the caller can drop the line entirely.
    """
    lam0 = float(line["lam0"])
    cls = get_profile_class(profile)
    if window is None and not cls.fitted:
        raise ValueError(f"profile {profile!r} needs an explicit window")
    if window is not None:
        fit_window, found = (float(window[0]), float(window[1])), True
    else:
        fit_window, found = auto_window(ref_spec, lam0, fallback=fallback,
                                        resolution=resolution)
    if not found:
        # no fittable core resolved: there is no trustworthy window, so don't measure an EW
        # at all — return it flagged (no Profile / NaN windows), to keep but exclude.
        nan = float("nan")
        return None, (nan, nan), (nan, nan), ("no_window",)
    sub = _crop(ref_spec, fit_window)
    xx, yy = sub.values
    try:
        prof = cls.measure(np.asarray(xx, float), np.asarray(yy, float), x0=lam0,
                           xunit=sub.xunit)
    except ValueError as exc:
        raise RuntimeError(f"profile fit failed at {lam0:.3f}") from exc
    if not cls.fitted:
        # an integrated width: the window is the measurement, nothing to centre or flag
        return prof, fit_window, fit_window, ()
    center = float(prof.params[0])
    hwhm = 0.5 * float(prof.fwhm)
    # synth window: SYNTH_HWHM_MULT HWHM about the fitted centre, but at least FIT_MARGIN_HWHM
    # HWHM beyond the fit window on each side so the fit region is always enclosed.
    if window is not None:
        synth_lo, synth_hi = fit_window
    else:
        synth_lo = min(center - SYNTH_HWHM_MULT * hwhm, fit_window[0] - FIT_MARGIN_HWHM * hwhm)
        synth_hi = max(center + SYNTH_HWHM_MULT * hwhm, fit_window[1] + FIT_MARGIN_HWHM * hwhm)
    flags = []
    if not (fit_window[0] <= center <= fit_window[1]):
        flags.append("off_center")
    return prof, (synth_lo, synth_hi), fit_window, tuple(flags)


def rew_gradient(synth_fn, a0, step, window=None, profile="gauss", x0=None,
                 batch_synth_fn=None):
    """Numerical ``dREW/dA`` at ``a0`` by central difference of the reduced EW over ±``step``.

    ``synth_fn(a) -> Spectrum`` is the full synthesis at abundance ``a`` (caller's convention).
    The reduced EW is measured over ``window`` (the fit window) with ``profile`` and reference
    ``x0`` using ``get_weq``'s default REW return — no manual log conversion. Returns the gradient
    per dex (``NaN`` if either flank cannot be measured).

    ``batch_synth_fn([a_hi, a_lo]) -> [Spectrum, Spectrum]`` (optional): synthesise both flanks
    in one call so they can be fanned across a single process pool; identical result to the two
    serial ``synth_fn`` calls, just scheduled together.
    """
    step = abs(float(step)) or 0.05
    try:
        if batch_synth_fn is not None:
            hi, lo = batch_synth_fn([a0 + step, a0 - step])
        else:
            hi, lo = synth_fn(a0 + step), synth_fn(a0 - step)
        rew_hi = _rew(hi, window, profile, x0)
        rew_lo = _rew(lo, window, profile, x0)
    except Exception:
        return float("nan")
    if not (np.isfinite(rew_hi) and np.isfinite(rew_lo)):
        return float("nan")
    return float((rew_hi - rew_lo) / (2.0 * step))


# sigma_a_from_rew moved to tumult.abundance (re-exported above).


def fit_line(line, *, synth_only_fn, synth_all_fn, synth_excl_fn, abund_fit_fn,
             ew_meas, a_guess=0.0, profile="gauss", weak_frac=0.10, blend_max=0.30,
             grad_step=0.125, grad_batch_fn=None, synth_full_fn=None):
    """Stage 2 (direct fit): screen the line, fit its abundance + vshift/vdop, and get ``dREW/dA``.

    Injected callables (synthesis-agnostic so this stays engine-independent and testable):

    * ``synth_only_fn(a) -> Spectrum`` — the line at abundance ``a`` with *only* the target
      species (``only_elements``); its fit-window width EW is the weak-line screen.
    * ``synth_all_fn(a) -> Spectrum`` — the synthesis the abundance is fitted with at abundance
      ``a`` (all species, or the target species alone when background lines are switched off),
      used for the ``dREW/dA`` gradient.
    * ``synth_excl_fn() -> Spectrum`` — the window with the target species removed
      (``exclude_elements``); with ``synth_full_fn(a_guess)`` gives the blend fraction.
    * ``synth_full_fn(a) -> Spectrum`` — the synthesis with every species, for the blend metric.
      The observation carries its blends whatever the fit synthesises. Defaults to
      ``synth_all_fn``.
    * ``abund_fit_fn(guess) -> (best: dict, fitted: Spectrum)`` — the adaptive fit of the abundance
      axis plus vshift/vdop against the observation over the fit window. ``best`` must expose the
      abundance under key ``"a"`` (plus ``"vshift"``/``"vdop"`` when fitted).

    ``ew_meas`` is the measured width EW (mÅ) from stage 1; ``line.fit_window`` bounds every
    screen/gradient EW. Mutates ``line`` in place and returns it. A screen rejection sets a flag,
    ``include=False`` and returns early (no fit).
    """
    fw, x0 = line.fit_window, line.lam0

    # 1) weak-line screen: species-only width EW in the fit window vs the measured EW (a ratio).
    ew_species = _ew(synth_only_fn(a_guess), fw, profile, x0)
    line.ew_species = float(ew_species)
    if not np.isfinite(ew_species) or ew_species < weak_frac * ew_meas:
        line.flags = ("weak_line",)
        line.include = False
        return line

    # 2) blend screen: contaminant width fraction in the fit window.
    blend, _, _ = blend_fraction((synth_full_fn or synth_all_fn)(a_guess), synth_excl_fn(), window=fw,
                                 profile=profile, x0=x0)
    line.blend_frac = float(blend)
    if np.isfinite(blend) and blend > blend_max:
        line.flags = ("blended",)
        line.include = False
        return line

    # 3) direct fit of abundance + vshift + vdop against the observation over the fit window.
    best, fitted = abund_fit_fn(a_guess)
    line.a = float(best["a"])
    line.vshift = float(best.get("vshift", np.nan))
    line.vdop = float(best.get("vdop", np.nan))
    line.fit_spec = fitted
    line.method = "fit"

    # 4) numerical dREW/dA around the best-fit abundance (full synthesis, fit window). Both
    # flanks are independent, so they go through the batch path (one shared pool) when given.
    line.drewda = rew_gradient(synth_all_fn, line.a, grad_step, window=fw, profile=profile,
                               x0=x0, batch_synth_fn=grad_batch_fn)

    # 5) flag an implausibly large velocity shift: compare |vshift| against the line's HWHM,
    # converted to a velocity (HWHM_λ / λ0 · c). A best fit that had to slide the synthetic by
    # more than half a line width is suspect, so mark it and exclude it by default.
    line.flags = ()
    hwhm_kms = velocity_hwhm(line.fwhm, line.lam0)
    if np.isfinite(line.vshift) and np.isfinite(hwhm_kms) and abs(line.vshift) > hwhm_kms:
        line.flags = ("vshift",)
        line.include = False
    return line


def finalize_line(line, ref_spec, snr, cont_pct, profile="gauss"):
    """Stage 3: re-measure the EW with SNR/continuum errors and propagate to ``sigma_a``.

    Keeps ``line.a`` from the direct fit and converts the EW error to an abundance error through
    the stored REW gradient ``line.drewda`` (dREW/dA). No synthesis — fast to rerun whenever the
    SNR / continuum-% inputs change. Mutates ``line`` in place.
    """
    prof = measure_ew(
        ref_spec, window=line.fit_window, snr=snr, cont_pct=cont_pct, profile=profile, x0=line.lam0)
    line.ew, line.ew_err = prof.maa, prof.width_err_maa
    if not np.isfinite(getattr(line, "drewda", float("nan"))):
        return line
    sa, sap, sac = sigma_a_from_rew(prof, line.drewda)
    line.sigma_a = float(sa)
    line.sigma_a_phot = float(sap)
    line.sigma_a_cont = float(sac)
    return line


# --- the three-stage element workflow ------------------------------------------------------
# Line selection, the per-line synthesis, and the error budget, as three calls. They are the
# same three the GUI's abundance wizard makes, with the rendering left to the caller: what
# comes back is the :class:`AbundanceResult`, and ``progress`` is how a caller watches it
# fill in line by line during the long middle stage.

def measure_lines(mset, ref_spec, *, linelist, profile, snr, cont_pct):
    """Stage 1: measure a mask set's switched-on windows on the reference observation.

    The windows are taken exactly as curated — nothing is re-derived, added or dropped, and
    the set itself is not written to. Each window's dominant transition is identified from
    the line list, which is what supplies the atomic identity the later stages read.

    Raises :class:`ValueError` when the set holds nothing measurable.
    """
    import uuid

    from . import linediag
    from .linelist import load_master, strip_ion_stages

    try:
        master = load_master(linelist)
        catalog = master[strip_ion_stages(master.ion) == mset.element]
    except Exception:
        # the catalog only names the dominant transition of each window; it never touches
        # the mask or the measurement
        logging.exception("Could not read the line list for %s", mset.element)
        catalog = None
    lines = linediag.measure_masks(ref_spec, [m for m in mset.masks if m.include],
                                   profile=profile, snr=snr, cont_pct=cont_pct,
                                   catalog=catalog)
    for ln in lines:
        if not ln.ion:
            ln.ion = mset.element
    if not lines:
        raise ValueError(f"No measurable windows in “{mset.name}”.")
    result = AbundanceResult(element=mset.element, symbol=mset.element, lines=lines,
                             summary=aggregate([]), stage="lines", profile=profile,
                             uid=uuid.uuid4().hex[:8], snr=snr, cont_pct=cont_pct,
                             maskset_uid=str(mset.uid), maskset_name=mset.name)
    from .measurement import mask_key
    result.maskset_key = mask_key(mset)
    result.name = f"{mset.element} · {len(lines)} lines"
    return result


def run_abundance(request, result, observation, *, progress=None, cancel_token=None):
    """Stage 2: screen and fit every included line of ``result``, in place.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest`; ``observation`` the spectrum the
    lines were measured on. Each line gets the syntheses :func:`fit_line` needs — the
    species in isolation, the full window, the window with the element removed, and the
    adaptive abundance fit — plus the dREW/dA flanks, which go in as one batch so they fold
    onto a single run's abundance-set axis. A request whose ``run_vals`` carry
    ``only_elements`` (background lines switched off) fits and differentiates the element's
    lines alone; the blend screen still compares against every species. ``progress(result)`` is called as each line
    lands, so a caller can stream partial results; it must not raise.
    """
    from dataclasses import asdict
    result.synthesis = asdict(request)
    from .synthesis import workflows as synth
    from .synthesis.cancel import CancelledError

    report = progress or (lambda _result: None)
    prefs, grid = request.prefs, request.grid
    target, atmos_file = request.target, request.atmos_file
    base_rv = dict(request.run_vals)
    feh = request.feh
    relative = bool(base_rv.get("relative_abundances", False))
    solar = synth.abund_table(prefs)
    base_abund = str(base_rv.get("abundances", "")).strip()
    symbol = result.symbol
    from .synthesis.request import Snapshot3D
    common_result = isinstance(request.selection, Snapshot3D) or request.departure_field is not None

    def synth_window(rv_over, window):
        # an override of None drops that key from the base run values
        rv = dict(base_rv)
        rv.update(rv_over)
        rv = {k: v for k, v in rv.items() if k not in rv_over or rv_over[k] is not None}
        if common_result:
            return request.with_run_vals(**rv).with_masks(
                [{'fit': window, 'synth': window}]).execute(
                    cancel_token=cancel_token).spectrum
        return synth.compute(target, prefs, rv, *padded_window(window), grid=grid,
                             atmos_file=atmos_file, cancel_token=cancel_token)

    def abund_terms(a):
        return abundance_override(symbol, a, base=base_abund, relative=relative,
                                  feh=feh, solar=solar)

    for ln in result.lines:
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        if not ln.include:
            continue          # skip step-1-rejected / user-excluded lines (keeps flags)

        def synth_only_fn(a, _w=ln.window):
            return synth_window({"abundances": abund_terms(a),
                                 "relative_abundances": relative,
                                 "only_elements": symbol}, _w)

        def synth_all_fn(a, _w=ln.window):
            return synth_window({"abundances": abund_terms(a),
                                 "relative_abundances": relative}, _w)

        def synth_all_batch(a_list, _w=ln.window):
            rvs = [{**base_rv, "abundances": abund_terms(a),
                    "relative_abundances": relative} for a in a_list]
            if common_result:
                return [request.with_run_vals(**rv).with_masks(
                    [{'fit': _w, 'synth': _w}]).execute(cancel_token=cancel_token).spectrum
                        for rv in rvs]
            return synth.compute_batch(target, prefs, rvs, *padded_window(_w), grid=grid,
                                       atmos_file=atmos_file, cancel_token=cancel_token)

        def synth_full_fn(a, _w=ln.window):
            return synth_window({"abundances": abund_terms(a),
                                 "relative_abundances": relative,
                                 "only_elements": None}, _w)

        def synth_excl_fn(_w=ln.window):
            return synth_window({"exclude_elements": symbol, "only_elements": None}, _w)

        def abund_fit_fn(guess, _w=ln.window, _fw=ln.fit_window):
            if common_result:
                return synth.fit_request_line(observation, request, symbol, guess, _w, _fw,
                                              cancel_token=cancel_token)
            return synth.fit_line(observation, target, prefs, base_rv, symbol, guess,
                                  *padded_window(_w), _fw,
                                  grid=grid, atmos_file=atmos_file,
                                  cancel_token=cancel_token)

        try:
            fit_line(ln, synth_only_fn=synth_only_fn, synth_all_fn=synth_all_fn,
                     synth_excl_fn=synth_excl_fn, abund_fit_fn=abund_fit_fn,
                     ew_meas=ln.ew, a_guess=0.0, profile=result.profile,
                     grad_batch_fn=synth_all_batch, synth_full_fn=synth_full_fn)
        except CancelledError:
            raise
        except Exception:
            logging.exception("Abundance: fit for line %.3f failed", ln.lam0)
            ln.flags = ("synth_failed",)
        result.relative = relative
        result.teff, result.logg, result.feh, result.vmic = target[0], target[1], feh, target[3]
        report(result)
    return result


def finalize_result(result, ref_spec, snr, cont_pct):
    """Stage 3: fold the SNR and continuum errors into every fitted line, then aggregate.

    Only lines that were screened and fitted in stage 2 carry a REW gradient, and only those
    can have an abundance error; the rest keep what they have. No synthesis, so this is cheap
    to rerun whenever the SNR or continuum-% inputs change.
    """
    for ln in result.lines:
        if np.isfinite(ln.drewda):
            finalize_line(ln, ref_spec, snr, cont_pct, profile=result.profile)
    result.snr, result.cont_pct = snr, cont_pct
    result.summary = aggregate(result.lines)
    result.stage = "uncert"
    result.termination = "completed"
    return result


# --- NLTE abundance corrections against a 1D LTE reference ---------------------------------
# The counterpart to measure_lines / run_abundance for the case where the "observation" is a
# synthetic 1D LTE reference computed straight from a model atom. Stage 1 measures the atom's
# own lines on that reference; stage 2 builds an NLTE curve of growth per line and inverts it
# to the abundance reproducing the reference EW. The reference sits at [X/Fe]=0, so that
# abundance *is* the correction dA = A(NLTE) - A(LTE).

CORR_LIMIT = 1.0       # the curve of growth spans +/-CORR_LIMIT dex in [X/Fe] about the reference
CORR_MAX_ROUNDS = 8    # bracketing rounds before giving up (each round is a batch of NLTE solves)
MIN_FIT_POINTS = 7     # a fit region sampled by fewer reference points than this is too sparse
                       # to trust - flag "few_points" and exclude the line by default


@dataclass
class CorrectionOpts:
    """Thresholds and ladder shape for the NLTE-correction workflow.

    The default profile is the trapezoidal width over each line's window: with the windows
    coming from the model atom (``line_grid``) the measurement is the same integral on the
    reference and on every ladder rung, which a fitted shape on a blended or damped line is
    not. It does couple lines whose windows overlap. A fitted profile needs ``line_grid``
    off, or accepts the atom windows as its fit region."""
    profile: str = "trapz"
    ew_min: float = 5.0                      # maa - reference lines weaker than this are dropped
    limit: float = CORR_LIMIT               # dex; the rungs stay within +-limit
    max_rounds: int = CORR_MAX_ROUNDS
    min_fit_points: int = MIN_FIT_POINTS
    # generate the bound-bound grid and line windows from the atom's opacity profiles on
    # the reference model (tumult.line_quadrature) instead of the atom's own q0/qmax; the
    # windows then also set the line masks
    line_grid: bool = True
    # synthesise the configured line list (blends from the other species) together with
    # the atom, in the reference and in every NLTE rung; off, the atom alone is synthesised
    linelist: bool = False


def correction_summary(lines):
    """The element mean of a set of corrections.

    Each correction carries a formal error after stage 2 - the reference EW's fit-covariance
    error (model mismatch, nonzero even on a noiseless synthetic spectrum) propagated through
    the CoG slope - so the inverse-variance weighted :func:`aggregate` is the right mean.
    :func:`aggregate_unweighted` remains the fallback for results without per-line errors
    (stage 1 only, or restored from a pre-error session).
    """
    weighted = aggregate(lines)
    if weighted["n"] > 0:
        return weighted
    return aggregate_unweighted(lines)


def species_symbol(ion):
    """``(species, symbol)`` from a line-list ``ion`` label like ``"Fe I"``.

    The species keeps its ionisation stage, with whitespace normalised; the symbol is that
    species without the stage (:func:`tumult.linelist.strip_ion_stage`) — the element for an
    atomic species, which is what the correction workflow's abundance overrides need.
    """
    from .linelist import strip_ion_stage

    species = " ".join(str(ion).split())
    return species, strip_ion_stage(species)


def measure_reference_lines(request, atom, *, ref_grid, opts=None, cancel_token=None):
    """Stage 1: synthesise the 1D LTE reference from ``atom`` and measure its lines.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest` naming the target parameters and
    preferences; ``ref_grid`` the grid flagged as the 1D LTE reference, which is interpolated
    at those parameters. The model atom alone defines the spectrum and its wavelength
    coverage, so no window is imposed - the reference spans the atom's full synthesised range,
    and that run's ``linelist_out.txt`` supplies the candidate lines.

    With ``opts.line_grid`` (the default) the bound-bound wavelength grid and the line
    windows come from :func:`tumult.line_quadrature.prepare_line_quadrature` on the
    reference model nearest the target, and the reference is synthesised on that grid with
    the atom copy carrying the windows; each line's window is then its fit and synthesis
    window (its line mask). Where the windows of two atom lines overlap, and the reference
    holds nothing but the atom (``opts.linelist`` off), each such line's window is narrowed
    to the reference flux maxima nearest its centre (:func:`local_max_window`), so the EW
    integrates the line's own profile and not its neighbour's core; with the line list in,
    the flux maxima belong to other species as much as to the atom, so the opacity windows
    are kept as they are. Stage 2 reuses the same preparation, recorded in
    ``result.quadrature``, so LTE reference and NLTE ladder share one grid and one set of
    windows, and measures the rungs over the same (narrowed) fit windows.

    Returns an :class:`AbundanceResult` at stage ``"lines"``, carrying the reference spectrum
    as ``ref_spec``. Raises :class:`ValueError` when the reference yields nothing measurable.
    """
    import uuid

    from .synthesis import workflows as synth

    opts = opts or CorrectionOpts()
    profile, ew_min = str(opts.profile), float(opts.ew_min)
    prefs, target = request.prefs, request.target

    quadrature, windows, run_atom = None, {}, atom
    if opts.line_grid:
        quadrature = _prepare_quadrature(request, atom, ref_grid)
        windows = _mask_windows(quadrature)
        run_atom = quadrature.atom_file
    ref_spec, run_dirs = synth.compute_reference(
        target, prefs, run_atom, ref_grid, cancel_token=cancel_token, return_runs=True,
        line_grid=quadrature.grid_file if quadrature else None,
        with_linelist=bool(opts.linelist))
    llist = synth.read_reference_linelist(run_dirs)
    if llist is None or len(llist) == 0:
        raise ValueError("No linelist_out.txt was produced by the reference run.")
    # with the line list included the dump holds every species; the candidates are the
    # atom's own lines
    element = atom_element(atom)
    species = llist["ion"].map(lambda ion: species_symbol(ion)[1])
    cand = llist[species == element]
    if {"i", "j"} <= set(cand.columns):
        # a transition split into isotopic / hyperfine components by the run's isotope data
        # is one line: keep its strongest record, centred on the atom's own wavelength when
        # the atom copy is at hand
        if "strength" in cand.columns:
            cand = cand.sort_values("strength", ascending=False)
        cand = cand.drop_duplicates(subset=["i", "j"])
        if windows:
            centres = _mask_centres(quadrature)
            cand = cand.assign(lam0=[centres.get((int(i), int(j)), l0) for i, j, l0
                                     in zip(cand["i"], cand["j"], cand["lam0"])])
    cand = cand.sort_values("lam0")
    if cand.empty:
        raise ValueError("The reference run produced no atom lines.")
    if windows and not {"i", "j"} <= set(cand.columns):
        raise ValueError("The reference linelist_out.txt carries no level indices.")

    blended = overlapping_windows(windows) if windows and not opts.linelist else set()
    lines = []
    for row in cand.itertuples():
        species, _symbol = species_symbol(getattr(row, "ion", ""))
        window = windows.get((int(row.i), int(row.j))) if windows else None
        if windows and window is None:
            continue
        fit = window
        if (int(row.i), int(row.j)) in blended:
            fit = local_max_window(ref_spec, float(row.lam0), window)
        try:
            lprof, win, fit_window, flags = measure_line(
                ref_spec, {"lam0": float(row.lam0)}, profile=profile, window=fit)
        except Exception:
            continue
        if window is not None:
            win = window          # the synth window stays the atom's own (its qmax span)
        if "no_window" in flags or lprof is None or not np.isfinite(lprof.maa):
            continue
        ew, func = lprof.maa, lprof.func
        fwhm = float(lprof.fwhm) if lprof.fitted else float("nan")
        # below the EW threshold a line is no candidate at all, whatever else is flagged;
        # listing it excluded would only pad the table with the atom's negligible lines
        if ew < ew_min:
            continue
        # a fit region resolved from too few reference samples can't be fitted reliably;
        # flag it (and exclude it) so it doesn't skew the corrections.
        flags = list(flags)
        n_pts = len(ref_spec.crop(xmin=fit_window[0], xmax=fit_window[1]).values[0])
        if n_pts < int(opts.min_fit_points):
            flags.append("few_points")
        # keep the fit-covariance EW error: even on the noiseless synthetic reference the
        # residual-scaled covariance measures the profile/model mismatch, and it propagates
        # to a per-line sigma(dA) through the CoG slope in stage 2 - the formal error the
        # weighted element mean uses.
        ew_err = float(lprof.width_err_maa)
        # flag (and exclude) a line whose EW is more uncertain than the EW itself: the
        # fit-covariance error is the only EW uncertainty on the synthetic reference, so
        # ew_err > ew means the measurement is not trustworthy.
        if np.isfinite(ew_err) and ew > 0 and ew_err > ew:
            flags.append("ew_unc")
        lines.append(LineResult(
            lam0=float(row.lam0), ion=species, e_low=float(getattr(row, "evi", np.nan)),
            loggf=float(getattr(row, "loggf", np.nan)), ew=ew,
            ew_err=ew_err, fwhm=fwhm,
            mask=LineMask(lam0=float(row.lam0), lmin=fit_window[0], lmax=fit_window[1],
                          smin=win[0], smax=win[1], note=species),
            fit_func=func, flags=tuple(flags), include=not flags))
    if not lines:
        raise ValueError(f"No atom lines above {ew_min:.0f} mÅ in the reference.")

    symbol = species_symbol(lines[0].ion)[1]
    result = AbundanceResult(
        element=symbol, symbol=symbol, lines=lines, summary=correction_summary([]),
        stage="lines", profile=profile, uid=uuid.uuid4().hex[:8], relative=True)
    result.teff, result.logg, result.feh, result.vmic = target
    result.ref_spec = ref_spec
    result.atom_path = str(atom)
    result.quadrature = quadrature.to_dict() if quadrature else {}
    result.with_linelist = bool(opts.linelist)
    if opts.linelist:
        result.linelist_dump = str(Path(run_dirs[0]) / "linelist_out.txt")
    result.name = Path(atom).name
    return result


def atom_element(atom):
    """The element symbol a classic model atom declares on its first entry line."""
    for line in Path(atom).read_text().splitlines():
        text = line.strip()
        if text and not text.startswith("*"):
            return text.split()[0].capitalize()
    raise ValueError(f"{atom} declares no element")


def _prepare_quadrature(request, atom, ref_grid):
    """The bound-bound grid and atom copy for ``atom`` on the reference model nearest the
    request's target, generated once and cached (see :mod:`tumult.line_quadrature`)."""
    from .line_quadrature import prepare_line_quadrature
    from .synthesis import workflows as synth

    prefs, target = request.prefs, request.target
    atmos = synth.nearest_marcs_atmosphere(ref_grid, target, prefs)
    vmic = float(target[3]) if len(target) > 3 and float(target[3]) > 0 else 1.0
    return prepare_line_quadrature(synth.preparation_root(prefs), atom, atmos,
                                   vmic=vmic, feh=float(request.feh))


def _mask_centres(quadrature):
    """``{(i, j): lam0}`` line centres of a preparation's atom copy, in air AA."""
    from .line_quadrature import line_windows
    from .physics import vac2air

    return {key: float(vac2air(lam0))
            for key, (lam0, _lo, _hi) in line_windows(quadrature.atom_file).items()}


def _mask_windows(quadrature):
    """``{(i, j): (lo, hi)}`` line windows of a preparation, in air AA as the masks are."""
    from .line_quadrature import line_windows
    from .physics import vac2air

    return {key: (float(vac2air(lo)), float(vac2air(hi)))
            for key, (_lam0, lo, hi) in line_windows(quadrature.atom_file).items()}


def _solution_spec(ln):
    """The line's rungs interpolated at its inverted abundance, over its synth window."""
    offsets = sorted(ln.cog_specs)
    return cog_spectrum(offsets, [ln.cog_specs[a] for a in offsets], ln.a)


def refit_line_window(result, ln, lo, hi, *, min_fit_points=MIN_FIT_POINTS):
    """Re-read one line after its fit window was moved to ``[lo, hi]`` (air Å), in place.

    The window is clamped to the line's synth window, which stays what the atom set. The
    reference EW is measured again over it; when the line carries its rungs
    (``cog_specs``, from :func:`run_nlte_corrections`) the curve of growth is re-read over
    the same window and inverted again, so the correction follows the edit without another
    solve. The reference-stage flags (``few_points``, ``ew_unc``) and the inversion's
    ``out_of_range`` are re-derived; the user's include toggle is kept unless the line
    falls out of range. The result's summary is refreshed.
    """
    smin, smax = ln.window
    if np.isfinite(smin) and np.isfinite(smax):
        lo, hi = max(float(lo), smin), min(float(hi), smax)
    if not hi > lo:
        return ln
    ln.mask.lmin, ln.mask.lmax = float(lo), float(hi)
    flags = []                      # every flag is re-derived below
    prof = measure_profile(result.ref_spec, (lo, hi), profile=result.profile, x0=ln.lam0)
    ln.ew, ln.ew_err = float(prof.maa), float(prof.width_err_maa)
    ln.fwhm = float(prof.fwhm) if prof.fitted else float("nan")
    ln.fit_func = prof.func
    n_pts = len(_crop(result.ref_spec, (lo, hi)).values[0])
    if n_pts < int(min_fit_points):
        flags.append("few_points")
    if np.isfinite(ln.ew_err) and ln.ew > 0 and ln.ew_err > ln.ew:
        flags.append("ew_unc")
    if ln.cog_specs:
        renorm = bool(getattr(ln.mask, "renorm", False))
        offsets = sorted(ln.cog_specs)
        ews = [_ew(ln.cog_specs[a], (lo, hi), result.profile, ln.lam0, renorm)
               for a in offsets]
        ln.a_grid, ln.ew_grid = tuple(offsets), tuple(float(w) for w in ews)
        if not (min(ews) <= ln.ew <= max(ews)):
            ln.a, ln.sigma_a, ln.fit_spec = float("nan"), float("nan"), None
            flags.append("out_of_range")
            ln.include = False
        else:
            ew_err = ln.ew_err if np.isfinite(ln.ew_err) else 0.0
            inv = cog_invert(offsets, ews, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
            ln.a = float(np.clip(inv["a"], offsets[0], offsets[-1]))
            ln.sigma_a = float(inv["sigma_a"])
            ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
            ln.fit_spec = _solution_spec(ln)
    ln.flags = tuple(flags)
    result.summary = correction_summary(result.lines)
    return ln


def run_nlte_corrections(request, result, atom, *, opts=None, progress=None,
                         cancel_token=None):
    """Stage 2: invert an NLTE curve of growth to each line's correction, in place.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest` naming the *NLTE target*
    atmosphere; ``result`` the stage-1 :class:`AbundanceResult` whose lines carry the
    reference EWs measured at [X/Fe] = 0.

    Every rung of the curve of growth is a full NLTE solve at its own abundance offset, and
    all lines share the rungs: one run per offset synthesises the atom's whole range, and each
    line's EW is measured over its own window on it. The offsets lie on the preferences'
    abundance lattice (``abund_step``) within ``+-opts.limit``. Rounds: the initial set is
    {-limit, 0, +limit}; each line's guess is the PCHIP inversion of the rungs computed so
    far; whatever lattice points bracket a guess and are not yet computed are solved next,
    and the loop ends when every line's guess is bracketed by computed rungs (or after
    ``opts.max_rounds``). A line whose reference EW lies outside the CoG at ``+-limit`` is
    flagged ``out_of_range`` and requests nothing more.

    ``progress(result)`` is called after every round, so a caller can stream the current
    guesses; it must not raise. Returns ``result`` at stage ``"corr"`` with its summary
    aggregated.
    """
    from dataclasses import asdict
    result.synthesis = asdict(request)
    from .synthesis import workflows as synth
    from .synthesis.cancel import CancelledError
    from .synthesis.request import Snapshot3D

    opts = opts or CorrectionOpts()
    report = progress or (lambda _result: None)
    prefs, grid = request.prefs, request.grid
    target, atmos_file = request.target, request.atmos_file
    symbol = result.symbol
    step = float(prefs.get("abund_step") or 0.25)
    if not step > 0:
        raise ValueError("abund_step must be positive")
    common_result = isinstance(request.selection, Snapshot3D) or request.departure_field is not None

    # the grid and atom copy stage 1 prepared, so both stages see the same line windows
    run_atom, line_grid = atom, None
    if result.quadrature:
        from .line_quadrature import PreparedLineQuadrature
        quadrature = PreparedLineQuadrature.from_dict(result.quadrature)
        if not (Path(quadrature.grid_file).is_file() and Path(quadrature.atom_file).is_file()):
            raise FileNotFoundError("The reference stage's bound-bound grid is gone; "
                                    "recompute the reference.")
        run_atom, line_grid = quadrature.atom_file, quadrature.grid_file

    def _rv(a):
        rv = {"mode": "NLTE", "only_atom": not result.with_linelist, "atom_lines": True,
              "atoms": [run_atom], "abundances": f"{symbol}={a:.4f}",
              "relative_abundances": True}
        if line_grid:
            rv["line_grid"] = line_grid
        return rv

    # set_spectrum=False leaves spectrum_params without R/daa and no window is given, so
    # every rung is the atom's whole range on its own frequency grid (the reference's
    # sampling); fold=False keeps each rung a departure solve of its own.
    def synth_offsets(offsets):
        if common_result:
            return [request.with_run_vals(**_rv(a)).with_masks([]).execute(
                cancel_token=cancel_token).spectrum for a in offsets]
        return synth.compute_batch(target, prefs, [_rv(a) for a in offsets], grid=grid,
                                   atmos_file=atmos_file, set_spectrum=False,
                                   cancel_token=cancel_token, fold=False)

    def snap(x):
        return float(round(float(np.round(x / step) * step), 10))

    limit = snap(max(step, np.floor(float(opts.limit) / step + 1e-9) * step))
    rungs = {}                         # offset -> full-range Spectrum

    def solve(offsets):
        offsets = sorted(set(snap(a) for a in offsets) - set(rungs))
        if not offsets:
            return
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        for a, spec in zip(offsets, synth_offsets(offsets)):
            rungs[a] = spec

    def widths(ln, offsets):
        renorm = bool(getattr(ln.mask, "renorm", False))
        return np.array([_ew(rungs[a], ln.fit_window, result.profile, ln.lam0, renorm)
                         for a in offsets])

    lines = [ln for ln in result.lines if ln.include]
    for ln in lines:
        ln.method = "cog"
    solve([-limit, 0.0, limit])
    for _round in range(int(opts.max_rounds)):
        needed = set()
        offsets = sorted(rungs)
        for ln in lines:
            if not ln.include:
                continue
            try:
                ews = widths(ln, offsets)
                ln.a_grid, ln.ew_grid = tuple(offsets), tuple(float(w) for w in ews)
                lo_ew, hi_ew = ews[offsets.index(-limit)], ews[offsets.index(limit)]
                if not (min(lo_ew, hi_ew) <= ln.ew <= max(lo_ew, hi_ew)):
                    ln.a, ln.sigma_a = float("nan"), float("nan")
                    ln.flags, ln.include = ("out_of_range",), False
                    continue
                ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                inv = cog_invert(offsets, ews, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
                ln.a = float(np.clip(inv["a"], -limit, limit))
                ln.sigma_a = float(inv["sigma_a"])
                ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
            except CancelledError:
                raise
            except Exception:
                logging.exception("Correction: line %.3f failed", ln.lam0)
                ln.flags, ln.include = ("synth_failed",), False
                continue
            # the lattice cell holding the guess; a guess on a lattice point is bracketed
            # by that point alone
            k = ln.a / step
            lo = float(np.clip(np.floor(k + 1e-9) * step, -limit, limit))
            hi = float(np.clip(np.ceil(k - 1e-9) * step, -limit, limit))
            needed.update(a for a in (lo, hi) if snap(a) not in rungs)
        report(result)
        if not needed:
            break
        solve(needed)
    else:
        logging.warning("Correction: %d rounds without bracketing every line", opts.max_rounds)

    # every line keeps its rungs over its synth window, so a later fit-window edit can
    # re-read the curve of growth without another solve (:func:`refit_line_window`); the
    # overlay for the jump view is the two bracketing rungs interpolated at the solution
    for ln in lines:
        lo, hi = ln.window
        if np.isfinite(lo) and np.isfinite(hi):
            ln.cog_specs = {a: rungs[a].crop(xmin=lo, xmax=hi) for a in ln.a_grid}
        if not ln.include or not np.isfinite(ln.a):
            continue
        ln.fit_spec = _solution_spec(ln)

    result.stage = "corr"
    result.summary = correction_summary(result.lines)
    return result
