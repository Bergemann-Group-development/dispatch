"""Headless, uncertainty-aware abundance and curve-of-growth workflows.

The synthesis-agnostic numerics (``cog_invert``, ``aggregate``, ``reduced_ew``, ``profile_fwhm``)
now live in :mod:`tumult.abundance` and are re-exported here for back-compat. This module keeps the
synthesis-driven, application-facing pieces:

* :func:`measure_lines` — stage 1: measure a mask set's windows on the reference (EW, FWHM).
* :func:`fit` — stage 2: screen the lines (weak / blended), then determine their abundances on
  shared lattice rungs, each rung one synthesis over every window: by curve-of-growth inversion
  (:func:`adaptive_cog`) or by **directly fitting** the synthetic to the observation (abundance +
  vshift + vdop, :func:`fit_profiles`). Either way ``dREW/dA`` is read off the same rungs.
* :func:`finalize_line` / :func:`aggregate` — stage 3: fold SNR + continuum-% errors into each
  line's EW, convert to an abundance error through the REW gradient, and combine into an
  inverse-variance element mean.

Error model: SNR noise is independent per pixel and the continuum offset is independent per line, so
both fold into each line's EW error inside ``tumult.utils.get_weq`` (``sigma`` → photon term,
``syst`` → continuum term). The abundance error works in *reduced* EW (``log10 W/λ``, unit
independent): ``σ_REW = σ_W/(W ln10)`` and ``σ_A = σ_REW / |dREW/dA|``.
"""
from __future__ import annotations

from copy import deepcopy
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
    auto_fit_window: tuple = ()  # the fit window the line was computed over, for a reset
    renorm: bool = False        # local continuum belongs to this target measurement
    continuum: tuple = ()      # (wavelength AA, prepared-target flux) nodes
    cog_step: float = float("nan")
    cog_origin: float = 0.0     # origin of the abundance lattice, in a_grid's convention
    missing_rungs: tuple = ()   # nearest lattice rungs needed after an EW edit
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
    ssr: float = float("nan")   # profile fit: summed squared flux residuals at the solution
    n_fit: int = 0              # ... over this many pixels (-> reduced_chi2)
    method: str = "fit"         # "fit" (direct) or "cog"
    flags: tuple = ()
    include: bool = True        # user toggle: count this line in the element mean

    # --- window views: the fit/synth windows now live on the embedded LineMask, but the
    # pipeline (measure_ew, the abundance, balance and correction fits) still
    # reads them as ``ln.fit_window`` / ``ln.window`` tuples, so expose them read-only here.
    @property
    def fit_window(self) -> tuple:
        return self.mask.fit_window

    @property
    def window(self) -> tuple:
        return self.mask.synth_window

    # --- session persistence (see AbundanceResult / the session snapshot) -------------
    _TUPLE_FIELDS = ("a_grid", "ew_grid", "flags", "continuum", "missing_rungs",
                     "auto_fit_window")

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
        if d.get("continuum"):
            d["continuum"] = tuple(tuple(p) for p in d["continuum"])
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
    """Explicit automatic local normalisation for standalone EW/CoG callers.

    Measurement workflows instead use their fixed nodes through ``measurement_spectrum``.
    """
    if not renorm:
        return spec
    from .continuum import divide_local
    return divide_local(spec)


def line_synthesis(line, spec):
    """A synthetic spectrum as ``line`` compares it with its target.

    With the line's local continuum enabled, the synthesis may carry a pseudo-continuum of
    its own (the far wings of lines outside the window), so it is divided by its own
    straight continuum, fitted to its upper envelope over the synth window
    (:func:`~tumult.continuum.divide_local`), the estimator that seeds the target's nodes.
    Stored rungs stay raw, so switching the continuum off restores the plain comparison.
    """
    if not line.renorm:
        return spec
    from .continuum import divide_local
    return divide_local(spec.convert_units(xunit=u.AA), window=line.window)


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
            prof = measure_ew(measurement_spectrum(spec, ln), window=ln.fit_window,
                              snr=snr, cont_pct=cont_pct, profile=profile, x0=ln.lam0)
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
    match an explicitly automatically normalised reference,
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
    method: str = ""                   # "cog", "profile" or "" for a staged result
    abundances: dict = field(default_factory=dict)  # composition the fit held, as Abundances
    provenance: dict = field(default_factory=dict)  # atmospheres, departures, quadrature

    def to_abundances(self):
        """The fitted composition as :class:`~tumult.synthesis.inputs.Abundances`.

        The composition the fit held fixed, with this element replaced by its mean, so
        ``request.update(abundances=result)`` synthesises at the fitted abundance. An NLTE
        correction is an abundance difference and has no such form.
        """
        from .synthesis.inputs import Abundances
        if self.stage == "corr":
            raise ValueError("NLTE corrections are abundance differences, not abundances.")
        mean = float(self.summary.get("a_mean", float("nan")))
        if not np.isfinite(mean):
            raise ValueError(f"No {self.symbol} abundance was determined.")
        if self.abundances:
            base = Abundances.from_dict(self.abundances)
        else:
            base = Abundances({}, "[X/Fe]" if self.relative else "A(X)")
        return base.with_value(self.symbol, mean)

    def stale_reasons(self, target=None):
        """:meth:`MeasuredLines.stale_reasons`, plus the lines whose edited local continuum
        left them without an abundance (``refit_required`` or ``missing_rungs``): the mean
        omits them until they are refitted or their rungs computed. Including or excluding a
        line needs no entry, because the mean is aggregated again from the per-line values."""
        reasons = super().stale_reasons(target)
        pending = sum(1 for ln in self.lines
                      if "refit_required" in ln.flags or ln.missing_rungs)
        if pending:
            reasons["continuum"] = (f"edited local continuum: {pending} line"
                                    f"{'s' if pending != 1 else ''} without an abundance")
        return reasons


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


def finalize_line(line, ref_spec, snr, cont_pct, profile="gauss"):
    """Stage 3: re-measure the EW with SNR/continuum errors and propagate to ``sigma_a``.

    Keeps ``line.a`` from stage 2 and converts the EW error to an abundance error through
    the stored REW gradient ``line.drewda`` (dREW/dA). No synthesis — fast to rerun whenever the
    SNR / continuum-% inputs change. Mutates ``line`` in place.
    """
    prof = measure_ew(
        measurement_spectrum(ref_spec, line), window=line.fit_window, snr=snr,
        cont_pct=cont_pct, profile=profile, x0=line.lam0)
    line.ew, line.ew_err = prof.maa, prof.width_err_maa
    if not np.isfinite(getattr(line, "drewda", float("nan"))):
        return line
    sa, sap, sac = sigma_a_from_rew(prof, line.drewda)
    line.sigma_a = float(sa)
    line.sigma_a_phot = float(sap)
    line.sigma_a_cont = float(sac)
    return line


# --- the three-stage element workflow ------------------------------------------------------
# Line selection (measure_lines), the lattice synthesis (fit, with ``measured=``) and the error
# budget (finalize_result), as the GUI's abundance wizard runs them; what comes back is the
# :class:`AbundanceResult`, and ``progress`` is how a caller watches it fill in round by round.

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
# The counterpart to measure_lines / fit for the case where the "observation" is a
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


def measure_reference_lines(request, atom, *, opts=None, quadrature=None, cancel_token=None):
    """Stage 1: synthesise the 1D LTE reference from ``atom`` and measure its lines.

    ``request`` is the :class:`~tumult.synthesis.SynthRequest` of the reference: its 1D
    atmosphere (a grid is interpolated at its star), star and composition; its departures
    do not enter. The model atom alone defines the spectrum and its wavelength coverage, so
    no window is imposed - the reference spans the atom's full synthesised range, and that
    run's ``linelist_out.txt`` supplies the candidate lines.

    With ``opts.line_grid`` (the default) the bound-bound wavelength grid and the line
    windows come from :func:`tumult.line_quadrature.prepare_line_quadrature` on the
    reference's MARCS model nearest the star (or ``quadrature``, prepared by the caller), and the
    reference is synthesised on that grid with
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

    from .departures import DepartureRequest
    from .synthesis import workflows as synth

    opts = opts or CorrectionOpts()
    profile, ew_min = str(opts.profile), float(opts.ew_min)

    windows, run_atom = {}, atom
    if opts.line_grid and quadrature is None:
        quadrature, _atmosphere = DepartureRequest((atom,)).prepare(request)
    if not opts.line_grid:
        quadrature = None
    if quadrature is not None:
        windows = _mask_windows(quadrature)
        run_atom = quadrature.atom_file
    ref_spec, run_dirs = synth.compute_reference(
        request, run_atom, cancel_token=cancel_token, return_runs=True,
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
            auto_fit_window=(float(fit_window[0]), float(fit_window[1])),
            fit_func=func, flags=tuple(flags), include=not flags))
    if not lines:
        raise ValueError(f"No atom lines above {ew_min:.0f} mÅ in the reference.")

    symbol = species_symbol(lines[0].ion)[1]
    result = AbundanceResult(
        element=symbol, symbol=symbol, lines=lines, summary=correction_summary([]),
        stage="lines", profile=profile, uid=uuid.uuid4().hex[:8], relative=True)
    star = request.star
    result.teff, result.logg, result.feh = star.teff, star.logg, star.feh
    result.vmic = float("nan") if star.vmic is None else star.vmic
    result.request = request
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
    return cog_spectrum(offsets, [line_synthesis(ln, ln.cog_specs[a]) for a in offsets],
                        ln.a)


def refit_line_window(result, ln, lo, hi, *, min_fit_points=MIN_FIT_POINTS):
    """Re-read one line after its fit window was moved to ``[lo, hi]`` (air Å), in place.

    The window is clamped to the line's synth window, which stays what the atom set. The
    reference EW is measured again over it; when the line carries its rungs
    (``cog_specs``, from :func:`run_nlte_corrections`) the curve of growth is re-read over
    the same window and inverted again, so the correction follows the edit without another
    solve. The reference-stage flags (``few_points``, ``ew_unc``) and the inversion's
    ``out_of_range`` are re-derived. The user's include toggle is kept, except that a line
    falling out of range is excluded and one coming back into range after such an
    exclusion is included again. The result's summary is refreshed.
    """
    smin, smax = ln.window
    if np.isfinite(smin) and np.isfinite(smax):
        lo, hi = max(float(lo), smin), min(float(hi), smax)
    if not hi > lo:
        return ln
    ln.mask.lmin, ln.mask.lmax = float(lo), float(hi)
    was_out = "out_of_range" in ln.flags
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
        offsets = sorted(ln.cog_specs)
        ews = [_ew(line_synthesis(ln, ln.cog_specs[a]), (lo, hi), result.profile, ln.lam0)
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
            if was_out:
                ln.include = True
    ln.flags = tuple(flags)
    result.summary = correction_summary(result.lines)
    return ln


def window_edited(ln):
    """Whether ``ln``'s fit window differs from the one it was computed over."""
    return bool(ln.auto_fit_window) and tuple(ln.fit_window) != tuple(ln.auto_fit_window)


def reset_line_window(result, ln):
    """Move ``ln``'s fit window back to the one it was computed over, in place.

    :func:`refit_line_window` re-reads the stored reference and rungs over that window, so
    the line's EW, curve of growth and correction are those of the computation and no
    synthesis runs.
    """
    if not ln.auto_fit_window:
        raise ValueError("This line does not record the window it was computed over.")
    return refit_line_window(result, ln, *ln.auto_fit_window)


def element_xfe(abundances, element, *, relative, feh, solar):
    """The [X/Fe] an ``abundances`` override string gives ``element`` (0 when absent)."""
    prefix = f"{element}=".lower()
    for term in str(abundances or "").replace(",", " ").split():
        if term.lower().startswith(prefix):
            value = float(term.split("=", 1)[1])
            return value if relative else abund_to_xfe(value, element, feh, solar)
    return 0.0


def lattice_rounds(lines, synthesize, solve_line, *, step, limit, max_rounds,
                   cancel_token=None, report=None, label="Abundance lattice", also=None):
    """Bracket every included line's abundance by its nearest rungs on a shared lattice.

    ``synthesize(offsets)`` returns, for each abundance offset [dex], the synthetic spectra
    aligned with ``lines``; an offset is only ever requested once, and all lines share it.
    ``solve_line(line, offsets, spectra)`` sets the line's solution from the computed rungs
    (sorted ``offsets``, the line's own ``spectra`` aligned with them) and returns its
    offset, or ``None`` after it has flagged and excluded the line.

    The offsets lie on a lattice of ``step`` dex within ``+-limit``. The first round
    synthesises ``(-limit, 0, +limit)``; each later one synthesises the lattice points
    bracketing every line's current solution that are not computed yet. A line is done when
    its solution lies between its nearest computed rungs. A line whose solve raises is
    flagged ``synth_failed`` and excluded. After ``max_rounds`` the last solutions stand.

    ``also()`` (optional) is called after each round's solves and returns further offsets
    whose lattice cells must be bracketed too (e.g. the same lines at other stars).

    Returns ``(rungs, bracketed)``: ``{offset: spectra}`` and whether every line was
    bracketed.
    """
    from .synthesis.cancel import CancelledError

    report = report or (lambda: None)
    step = float(step)
    limit = lattice_limit(step, limit)

    def snap(x):
        return float(round(float(np.round(x / step) * step), 10))

    rungs = {}                         # offset -> spectra aligned with lines

    def cell(a):
        """The uncomputed rungs of the lattice cell holding ``a``; one on a lattice point
        is bracketed by that point alone."""
        k = a / step
        lo = float(np.clip(np.floor(k + 1e-9) * step, -limit, limit))
        hi = float(np.clip(np.ceil(k - 1e-9) * step, -limit, limit))
        return {x for x in (lo, hi) if snap(x) not in rungs}

    def compute(offsets):
        offsets = sorted(set(snap(a) for a in offsets) - set(rungs))
        if not offsets:
            return
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        for a, specs in zip(offsets, synthesize(offsets)):
            rungs[a] = list(specs)

    active = [i for i, ln in enumerate(lines) if ln.include]
    compute([-limit, 0.0, limit])
    bracketed = False
    for _round in range(int(max_rounds)):
        needed = set()
        offsets = sorted(rungs)
        for i in active:
            ln = lines[i]
            if not ln.include:
                continue
            try:
                a = solve_line(ln, offsets, [rungs[o][i] for o in offsets])
            except CancelledError:
                raise
            except Exception:
                logging.exception("%s: line %.3f failed", label, ln.lam0)
                ln.flags = tuple(dict.fromkeys(ln.flags + ("synth_failed",)))
                ln.include = False
                continue
            if a is None:
                continue
            needed.update(cell(a))
        if also is not None:
            for a in also():
                if np.isfinite(a):
                    needed.update(cell(a))
        report()
        if not needed:
            bracketed = True
            break
        compute(needed)
    else:
        logging.warning("%s: %d rounds without bracketing every line", label, max_rounds)
    return rungs, bracketed


def lattice_limit(step, limit):
    """``limit`` snapped down onto the lattice of ``step`` dex, at least one step."""
    step = float(step)
    if not step > 0:
        raise ValueError("The abundance lattice step must be positive.")
    return float(round(float(max(step, np.floor(float(limit) / step + 1e-9) * step)), 10))


def _exclude(line, flag):
    line.a, line.sigma_a = float("nan"), float("nan")
    line.ssr, line.n_fit = float("nan"), 0
    line.flags = tuple(dict.fromkeys(line.flags + (flag,)))
    line.include = False


def adaptive_cog(lines, synthesize, *, step, limit, max_rounds, profile, cancel_token=None,
                 report=None, also=None):
    """Bracket every line's abundance on a lattice of curve-of-growth rungs, in place.

    ``lines`` carry the EW to reproduce (``ew``, ``ew_err``) and their fit windows.
    ``synthesize(offsets)`` returns, for each abundance offset [dex], the synthetic spectra
    aligned with ``lines`` (see :func:`lattice_rounds`).

    Each round inverts every line's curve of growth through the computed rungs (PCHIP in
    reduced EW, :func:`cog_invert`). A line whose EW lies outside the curve at ``+-limit``
    is flagged ``out_of_range`` and excluded.

    Sets ``a`` (the offset), ``sigma_a``, ``drewda``, ``a_grid``/``ew_grid``, ``cog_specs``
    and ``fit_spec`` on each line. Returns ``(rungs, bracketed)``: the offsets computed and
    whether every line was bracketed. ``also`` is passed to :func:`lattice_rounds`.
    """
    active = [ln for ln in lines if ln.include]
    for ln in active:
        ln.method = "cog"
        ln.cog_step, ln.cog_origin = float(step), 0.0
    edge = lattice_limit(step, limit)
    measured = {}          # (line, offset) -> EW: a rung's spectrum never changes

    def invert(ln, offsets, specs):
        for a, spec in zip(offsets, specs):
            if (id(ln), a) not in measured:
                measured[id(ln), a] = _ew(line_synthesis(ln, spec), ln.fit_window, profile,
                                          ln.lam0)
        ews = np.array([measured[id(ln), a] for a in offsets])
        ln.a_grid, ln.ew_grid = tuple(offsets), tuple(float(w) for w in ews)
        lo_ew, hi_ew = ews[offsets.index(-edge)], ews[offsets.index(edge)]
        if not (min(lo_ew, hi_ew) <= ln.ew <= max(lo_ew, hi_ew)):
            _exclude(ln, "out_of_range")
            return None
        ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
        inv = cog_invert(offsets, ews, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
        ln.a = float(np.clip(inv["a"], -edge, edge))
        ln.sigma_a = float(inv["sigma_a"])
        ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
        slope = float(inv["dadrew"])
        ln.drewda = 1.0 / slope if np.isfinite(slope) and slope != 0 else float("nan")
        return ln.a

    rungs, bracketed = lattice_rounds(
        lines, synthesize, invert, step=step, limit=limit, max_rounds=max_rounds,
        cancel_token=cancel_token, report=report, label="Curve of growth", also=also)

    # every line keeps its rungs over its synth window, so a later fit-window edit can
    # re-read the curve of growth without another synthesis (:func:`refit_line_window`);
    # the overlay is the two bracketing rungs interpolated at the solution
    index = {id(ln): i for i, ln in enumerate(lines)}
    for ln in active:
        lo, hi = ln.window
        if np.isfinite(lo) and np.isfinite(hi) and ln.a_grid:
            ln.cog_specs = {a: rungs[a][index[id(ln)]].crop(xmin=lo, xmax=hi)
                            for a in ln.a_grid}
        if ln.include and np.isfinite(ln.a):
            ln.fit_spec = _solution_spec(ln)
    return sorted(rungs), bracketed


PROFILE_FIT_PARAMS = 3          # abundance, vshift, vdop — the dof a reduced chi-square loses


def _fit_crop(line, observed):
    """The local-continuum crop arguments the profile fit evaluates its model through, as
    :meth:`tumult.spectrum.Spectrum.fit` assembles them from what it is handed."""
    nodes = line.continuum if line.renorm else None
    if not nodes:
        return {"continuum": None}
    return {"continuum": nodes,
            "continuum_window": (observed.xdata.min(), observed.xdata.max())}


def _record_residuals(ln, grid, observed):
    """Keep the objective the profile fit minimised, at the solution the line now carries.

    The optimiser's own last evaluation sits at unrounded parameters, while ``a``, ``vshift``
    and ``vdop`` are stored rounded to the precision they are shown at. Re-evaluating the
    same model there — one more interpolation of rungs already in hand, no synthesis — is
    what makes the reported residual describe the solution that is stored, drawn on the
    spectrum and scanned by :func:`chi2_profile`. Sets ``ssr`` and ``n_fit``.
    """
    model = grid.crop(newx=observed.xdata, a=ln.a, vshift=ln.vshift, vdop=ln.vdop,
                      **_fit_crop(ln, observed))
    residual = observed.values[1] - model.values[1]
    ln.ssr, ln.n_fit = float(np.sum(residual ** 2)), int(residual.size)


def _chi2_scale(n_fit, snr):
    """What turns a summed squared residual over ``n_fit`` pixels into a reduced chi-square
    at a continuum noise level of ``1/snr``. NaN when there is no scale to apply."""
    dof = int(n_fit) - PROFILE_FIT_PARAMS
    snr = float(snr) if snr else 0.0
    if dof <= 0 or snr <= 0:
        return float("nan")
    return snr ** 2 / dof


def reduced_chi2(line, snr):
    """``line``'s profile fit per degree of freedom, against a noise level of ``1/snr``.

    The fit itself is unweighted (:func:`_record_residuals`), so the SNR enters only here:
    the same continuum noise that sets the photon EW error. NaN for a line that was not
    profile-fitted, has no usable SNR, or has no degrees of freedom left.
    """
    return (float(getattr(line, "ssr", float("nan")))
            * _chi2_scale(getattr(line, "n_fit", 0), snr))


CHI2_SCAN_POINTS = 61       # samples of the abundance axis in a chi-square scan


def chi2_profile(line, target, snr, *, points=CHI2_SCAN_POINTS):
    """``line``'s fit objective along the abundance axis, as reduced chi-square.

    Evaluates exactly the model the fit did (:func:`_fit_line_profile`): the line's cached
    rungs interpolated in abundance, shifted and broadened by the fitted ``vshift`` and
    ``vdop`` and multiplied by its local continuum, compared with ``target`` inside the fit
    window. Only the abundance moves, so this is the fit's own objective sliced through its
    solution — which is what makes the reported minimum the one the optimiser stopped at.
    No synthesis runs: the rungs are the ones already cached on the line.

    Returns ``(abundances, chi2)`` spanning the cached rungs, with the solution itself
    among the samples; empty when the line carries no profile fit to scan.
    """
    from .spectrum_grid import SpectrumGrid

    if len(line.cog_specs) < 2 or not (np.isfinite(line.vshift) and np.isfinite(line.vdop)):
        return (), ()
    offsets = sorted(line.cog_specs)
    grid = SpectrumGrid([line_synthesis(line, line.cog_specs[a]) for a in offsets],
                        a=list(offsets))
    observed = target.crop(xmin=line.fit_window[0], xmax=line.fit_window[1])
    flux = observed.values[1]
    scale = _chi2_scale(flux.size, snr)
    crop = {"vshift": float(line.vshift), "vdop": float(line.vdop),
            **_fit_crop(line, observed)}
    scan = np.linspace(offsets[0], offsets[-1], points)
    if np.isfinite(line.a) and offsets[0] <= line.a <= offsets[-1]:
        scan = np.unique(np.append(scan, float(line.a)))
    chi2 = [float(np.sum((flux - grid.crop(newx=observed.xdata, a=float(a),
                                           **crop).values[1]) ** 2) * scale)
            for a in scan]
    return tuple(float(a) for a in scan), tuple(chi2)


def ew_window(line, profile):
    """The wavelength span ``line``'s equivalent width covers, or None when it has no end.

    An integrating profile (``"trapz"``, ``"pchip"``) measures the flux inside the fit
    window and nothing beyond it, so that window is the area. A fitted analytic profile
    integrates its own wings without limit, and no finite span holds all of that area —
    a Voigt's damping wings least of all. None says so, and leaves a caller drawing the
    area to show as much of it as it has room for.
    """
    if not get_profile_class(profile).fitted or not np.all(np.isfinite(line.fit_window)):
        return line.fit_window
    return None


def _fit_line_profile(ln, offsets, specs, target, *, guess):
    """Fit ``ln``'s abundance, ``vshift`` and ``vdop`` to ``target`` inside its fit window,
    interpolating its rungs ``specs`` at ``offsets`` along the abundance axis. A fixed local
    continuum (``renorm``) multiplies the model. Sets ``a``, ``vshift``, ``vdop`` and the
    fit's own residuals (:func:`_record_residuals`); returns ``(grid, best)``."""
    from .spectrum_grid import SpectrumGrid

    grid = SpectrumGrid(specs, a=list(offsets))
    observed = target.crop(xmin=ln.fit_window[0], xmax=ln.fit_window[1])
    best, _fitted = grid.fit(
        observed, fit=["a", "vshift", "vdop"],
        p0={"a": guess, "vshift": 0.0, "vdop": 1.0},
        continuum=ln.continuum if ln.renorm else None,
        xtol={"a": 0.01, "vshift": 0.1, "vdop": 0.1})
    ln.a = float(best["a"])
    ln.vshift, ln.vdop = float(best["vshift"]), float(best["vdop"])
    _record_residuals(ln, grid, observed)
    return grid, best


def _profile_solution(ln, grid, best, offsets, specs, profile):
    """The accepted profile fit's overlay (over the synth window, normalised) and ``drewda``,
    the slope of the PCHIP interpolant through the rungs' reduced EWs at the solution."""
    from .utils import pchip_mono

    ln.fit_spec = grid.crop(**best)
    rews = [_rew(spec, ln.fit_window, profile, ln.lam0) for spec in specs]
    ln.drewda = float(pchip_mono(np.asarray(offsets, float),
                                 np.asarray(rews, float)).derivative()(ln.a))


def fit_profiles(lines, synthesize, target, *, step, limit, max_rounds, profile,
                 cancel_token=None, report=None):
    """Fit every line's abundance, velocity shift and broadening on shared lattice rungs.

    ``lines`` carry their fit and synth windows (and a fixed local continuum when
    ``renorm``); ``synthesize(offsets)`` is as in :func:`lattice_rounds` and ``target`` is
    the observation. Each round interpolates every line's computed rungs along the
    abundance axis (monotone PCHIP per pixel, :class:`~tumult.spectrum_grid.SpectrumGrid`)
    and fits abundance, ``vshift`` and ``vdop`` to the target inside the fit window. A
    solution at ``+-limit`` is flagged ``out_of_range`` and excluded.

    ``drewda`` is the slope of the PCHIP interpolant through the rungs' reduced EWs at
    the solution, so the gradient and the fit rest on the same lattice. A line whose
    ``|vshift|`` exceeds its measured half width is flagged ``vshift`` and excluded.

    Sets ``a`` (the offset), ``vshift``, ``vdop``, ``drewda`` and ``fit_spec`` on each line,
    and keeps its rungs over its synth window (``cog_specs``, ``cog_step``, ``cog_origin``)
    so an edited local continuum can be fitted again without another synthesis
    (:func:`refit_line_profile`). Returns ``(rungs, bracketed)`` as :func:`adaptive_cog`
    does.
    """
    active = [ln for ln in lines if ln.include]
    for ln in active:
        ln.method = "fit"
        ln.cog_step, ln.cog_origin = float(step), 0.0
    edge = lattice_limit(step, limit) - 1e-3
    previous = {}                      # each line starts from its last round's solution

    def fit(ln, offsets, specs):
        specs = [line_synthesis(ln, spec) for spec in specs]
        grid, best = _fit_line_profile(ln, offsets, specs, target,
                                       guess=previous.get(id(ln), 0.0))
        previous[id(ln)] = ln.a
        if abs(ln.a) >= edge:
            _exclude(ln, "out_of_range")
            return None
        _profile_solution(ln, grid, best, offsets, specs, profile)
        return ln.a

    rungs, bracketed = lattice_rounds(
        lines, synthesize, fit, step=step, limit=limit, max_rounds=max_rounds,
        cancel_token=cancel_token, report=report, label="Profile fit")
    index = {id(ln): i for i, ln in enumerate(lines)}
    for ln in active:
        lo, hi = ln.window
        if np.isfinite(lo) and np.isfinite(hi):
            ln.cog_specs = {a: specs[index[id(ln)]].crop(xmin=lo, xmax=hi)
                            for a, specs in rungs.items()}
    for ln in active:
        hwhm_kms = velocity_hwhm(ln.fwhm, ln.lam0)
        if (ln.include and np.isfinite(ln.vshift) and np.isfinite(hwhm_kms)
                and abs(ln.vshift) > hwhm_kms):
            ln.flags = tuple(dict.fromkeys(ln.flags + ("vshift",)))
            ln.include = False
    return sorted(rungs), bracketed


def run_nlte_corrections(request, result, atom, *, opts=None, progress=None,
                         cancel_token=None):
    """Stage 2: invert an NLTE curve of growth to each line's correction, in place.

    ``request`` is the :class:`~tumult.synthesis.SynthRequest` naming the *NLTE target*
    atmosphere, star and composition (its departures do not enter: every rung solves
    ``atom``'s own); ``result`` the stage-1 :class:`AbundanceResult` whose lines carry the
    reference EWs measured at that composition.

    Every rung of the curve of growth is a full NLTE solve at its own abundance offset,
    applied on top of the element's abundance in the request's composition (the standard
    composition, [X/Fe] = 0, without an override); every other override is kept. All lines
    share the rungs: one run per offset synthesises the atom's whole range, and each line's
    EW is measured over its own window on it. The rounds are those of
    :func:`adaptive_cog`, with the lattice step ``abund_step`` and ``opts.limit``, so each
    line's correction ``A(NLTE) - A(LTE)`` ends bracketed by its nearest computed rungs.

    ``progress(result)`` is called after every round, so a caller can stream the current
    guesses; it must not raise. Returns ``result`` at stage ``"corr"`` with its summary
    aggregated and ``termination`` ``"completed"`` or ``"round_limit"``.
    """
    from .synthesis import workflows as synth
    from .synthesis.api import NodeStore
    from .synthesis.inputs import Snapshot
    from .synthesis.operations import run_values as complete_values

    opts = opts or CorrectionOpts()
    report = progress or (lambda _result: None)
    result.request = request
    request = request.update(departures=None)
    symbol = result.symbol
    step = request.config.abundance_step
    snapshot = isinstance(request.atmosphere, Snapshot)

    # the grid and atom copy stage 1 prepared, so both stages see the same line windows
    run_atom, line_grid = atom, None
    if result.quadrature:
        from .line_quadrature import PreparedLineQuadrature
        quadrature = PreparedLineQuadrature.from_dict(result.quadrature)
        if not (Path(quadrature.grid_file).is_file() and Path(quadrature.atom_file).is_file()):
            raise FileNotFoundError("The reference stage's bound-bound grid is gone; "
                                    "recompute the reference.")
        run_atom, line_grid = quadrature.atom_file, quadrature.grid_file

    base_rv = complete_values(request)
    relative = bool(base_rv.get("relative_abundances", True))
    base_abund = str(base_rv.get("abundances", "") or "")
    solar = synth.abund_table(request.config.to_prefs()) if not relative else {}
    feh = float(request.star.feh)
    base_x = element_xfe(base_abund, symbol, relative=relative, feh=feh, solar=solar)

    def _rv(a):
        rv = {"mode": "NLTE", "only_atom": not result.with_linelist, "atom_lines": True,
              "atoms": [run_atom],
              "abundances": abundance_override(symbol, base_x + a, base=base_abund,
                                               relative=relative, feh=feh, solar=solar),
              "relative_abundances": relative}
        if line_grid:
            rv["line_grid"] = line_grid
        return rv

    lines = [ln for ln in result.lines if ln.include]
    # the windows the rungs are measured over are what a later window reset returns to
    for ln in result.lines:
        if not ln.auto_fit_window:
            ln.auto_fit_window = tuple(float(v) for v in ln.fit_window)

    def coordinate(a):
        """A rung's store coordinate: the element's value as :func:`abundance_override`
        writes it."""
        x = base_x + a
        return round(float(x if relative else xfe_to_abund(x, symbol, feh, solar)), 4)

    # set_spectrum=False leaves spectrum_params without R/daa and no window is given, so
    # every rung is the atom's whole range on its own frequency grid (the reference's
    # sampling); each rung is a run, and a departure solve, of its own.
    store = NodeStore()

    def synth_offsets(offsets):
        if snapshot:
            specs = [synth.execute_request(request.update(windows=()), _rv(a),
                                           cancel_token=cancel_token).spectrum
                     for a in offsets]
        else:
            specs = synth.compute_batch(request, [_rv(a) for a in offsets],
                                        set_spectrum=False, cancel_token=cancel_token,
                                        store=store,
                                        abundances=[coordinate(a) for a in offsets])
        return [[spec] * len(lines) for spec in specs]

    _rungs, bracketed = adaptive_cog(
        lines, synth_offsets, step=step, limit=opts.limit, max_rounds=opts.max_rounds,
        profile=result.profile, cancel_token=cancel_token, report=lambda: report(result))

    result.stage = "corr"
    result.termination = "completed" if bracketed else "round_limit"
    result.summary = correction_summary(result.lines)
    return result


# --- complete analysis operations ----------------------------------------------------------
# The staged functions above let a caller inspect and edit each stage. The operations below
# own the whole lifecycle; the GUI's tools reuse a curated measurement through ``measured``.

def _target_spectrum(target):
    """A :class:`Spectrum` from a spectrum or a file path, and a provenance record."""
    source = ""
    if isinstance(target, (str, Path)):
        from .spectrum_io import load_file
        source = str(Path(target).expanduser())
        target = load_file(source)
    if not isinstance(target, Spectrum):
        raise TypeError("target must be a Spectrum or the path of a spectrum file.")
    from .measured import target_digest
    xx = np.asarray(target.values[0], dtype=float)
    if xx.size < 2:
        raise ValueError("The target spectrum has fewer than two samples.")
    digest = target_digest(target)
    record = {"source": source, "points": int(xx.size),
              "range": [float(np.min(xx)), float(np.max(xx))],
              "unit": str(getattr(target, "xunit", "") or ""), "sha1": digest}
    return target, record


def _line_mask_set(line_mask):
    """A :class:`~tumult.linemask.LineMaskSet` from a set or a mask file path."""
    from .linemask import LineMaskSet, parse_linemask_file
    if isinstance(line_mask, (str, Path)):
        line_mask = parse_linemask_file(Path(line_mask).expanduser())
    if not isinstance(line_mask, LineMaskSet):
        raise TypeError("line_mask must be a LineMaskSet or the path of a mask file.")
    if not line_mask.element:
        raise ValueError("The line mask set names no element; set its `element`.")
    if not any(mask.include for mask in line_mask.masks):
        raise ValueError(f"“{line_mask.name}” has no included windows.")
    return line_mask


def _check_method(method):
    if method not in ("cog", "profile"):
        raise ValueError(f"Unknown abundance method {method!r}; use 'cog' (adaptive curve of "
                         "growth) or 'profile' (direct profile fit).")


class _LatticeContext:
    """The syntheses of one lattice analysis (curve of growth or profile fit) on one request.

    Resolves the request once: the model atom's quadrature for requested departures, the
    element's abundance in the composition and the windows every line is synthesised over.
    ``variant`` builds the run values of one synthesis; every synthesis covers all lines'
    windows in one run through a window file.

    ``store`` (a :class:`~tumult.synthesis.api.NodeStore`) holds the rung runs of the
    whole workflow, so contexts at other stars reuse the atmosphere nodes they share; a
    context gets its own when none is given. The rungs are pruned per node at
    ``limit`` above the composition's abundance, and the first context on a store fixes
    that for the workflow.

    ``probes`` are further stellar parameters (mappings) at which every ``"fit"`` rung is
    also interpolated, from the same node runs and any further nodes the probes need;
    :func:`fit_lines_cog` re-inverts each line at them (``probe_solutions``), which is
    what a Jacobian in the stellar parameters is built from.
    """

    def __init__(self, request, lines, element, *, background_lines, limit, cancel_token,
                 store=None, probes=()):
        from .departures import DepartureGrid, Departures
        from .synthesis import workflows as synth
        from .synthesis.api import NodeStore
        from .synthesis.inputs import Snapshot
        from .synthesis.operations import departure_run_values, run_values

        self.departure_values, self.quadrature, self.quadrature_atmosphere = \
            departure_run_values(request)
        self.request, self.element = request, element
        self.solved_departures = isinstance(request.departures, (Departures, DepartureGrid))
        # snapshots and solved departures run through execute_request, which reads them;
        # the batch path schedules every variant's runs together on 1D atmospheres
        self.direct = self.solved_departures or isinstance(request.atmosphere, Snapshot)
        self.lte_request = request.update(departures=None)
        base = run_values(request, **self.departure_values)
        self.prefs = request.config.to_prefs()
        self.solar = synth.abund_table(self.prefs)
        self.feh = float(request.star.feh)
        self.relative = bool(base.get("relative_abundances", True))
        self.base_abund = str(base.get("abundances", "") or "")
        self.base_x = element_xfe(self.base_abund, element, relative=self.relative,
                                  feh=self.feh, solar=self.solar)
        self.background_lines = bool(background_lines)
        self.cancel_token = cancel_token
        self.lines = [ln for ln in lines
                      if np.isfinite(ln.window[0]) and np.isfinite(ln.window[1])]
        self.windows = [padded_window(ln.window) for ln in self.lines]
        self.prune = {"abundances": self._abundances(self.base_x + float(limit)),
                      "relative_abundances": self.relative}
        self.store = NodeStore() if store is None else store
        self.probes = [dict(probe) for probe in probes]
        self.probe_rungs = {}        # offset -> [probe][window] spectra
        self.probe_solutions = []    # [probe] {id(line): ([X/Fe], sigma)}
        self.probe_ews = {}          # (probe, id(line), offset) -> EW, measured once

    def _abundances(self, xfe):
        return abundance_override(self.element, xfe, base=self.base_abund,
                                  relative=self.relative, feh=self.feh, solar=self.solar)

    def _coordinate(self, xfe):
        """The element's value as :meth:`_abundances` writes it: the rung's store coordinate."""
        value = xfe if self.relative else xfe_to_abund(xfe, self.element, self.feh, self.solar)
        return round(float(value), 4)

    def variant(self, xfe, *, species="fit"):
        """Run values at [X/Fe] ``xfe``. ``species`` is ``"fit"`` (what the rungs
        synthesise), ``"all"``, ``"only"`` (the element alone) or ``"excluded"``."""
        # the element's term stays in every variant, so a single-variant batch always
        # folds its windows into one run, the excluded arm included
        values = {"abundances": self._abundances(xfe),
                  "relative_abundances": self.relative}
        if species == "fit":
            species = "all" if self.background_lines else "only"
        if species == "only":
            values["only_elements"] = self.element
        elif species == "excluded":
            values["exclude_elements"] = self.element
        return values

    def synthesize(self, variants, *, lte=False, rungs=None, probes=()):
        """``spectra[variant][window]`` for run-value ``variants``: LTE syntheses at the
        request's composition with ``lte``, else with the request's departures. ``rungs``
        gives the [X/Fe] of each ``"fit"`` variant, whose runs then go through the store.
        With ``probes`` the result is ``spectra[star][variant][window]``, the request's
        star first."""
        from .synthesis import workflows as synth
        if lte:
            request, extra = self.lte_request, {}
        else:
            request, extra = self.request, self.departure_values
        if self.direct:
            if probes:
                raise ValueError("Probe stars need a 1D atmosphere grid without solved "
                                 "departures.")
            windowed = request.update(windows=self.windows)
            return [list(synth.execute_request(windowed, {**extra, **variant},
                                               cancel_token=self.cancel_token).spectra)
                    for variant in variants]
        stored = {} if rungs is None else {
            "store": self.store, "abundances": [self._coordinate(x) for x in rungs]}
        return synth.compute_window_batch(
            request, [{**extra, **variant} for variant in variants], self.windows,
            prune=self.prune, cancel_token=self.cancel_token, probes=probes, **stored)

    def synthesize_offsets(self, lines):
        """``synthesize(offsets)`` for :func:`lattice_rounds`: the ``"fit"`` variant at each
        [X/Fe] offset from the composition, its spectra aligned with ``lines``."""
        index = {id(ln): w for w, ln in enumerate(self.lines)}
        x0 = self.base_x

        def synthesize(offsets):
            rungs = [x0 + a for a in offsets]
            spectra = self.synthesize([self.variant(x) for x in rungs], rungs=rungs,
                                      probes=self.probes)
            if self.probes:
                spectra, at_probes = spectra[0], spectra[1:]
                for vi, a in enumerate(offsets):
                    self.probe_rungs[a] = [probe[vi] for probe in at_probes]
            return [[specs[index[id(ln)]] for ln in lines] for specs in spectra]
        return synthesize

    def solve_probes(self, lines, profile):
        """Invert each solved line of ``lines`` through its curve of growth at every probe.

        Uses the rungs each line was bracketed with (``a_grid`` as offsets), so it runs
        after :func:`adaptive_cog` and before the offsets are shifted. Each solution is
        ``([X/Fe], sigma_a)`` as :func:`adaptive_cog` would set them; a line whose EW leaves
        a probe's curve gets NaN there.
        """
        index = {id(ln): w for w, ln in enumerate(self.lines)}
        self.probe_solutions = []
        for k in range(len(self.probes)):
            solutions = {}
            for ln in lines:
                if not (ln.include and np.isfinite(ln.a) and ln.a_grid):
                    continue
                offsets = list(ln.a_grid)
                for a in offsets:
                    if (k, id(ln), a) not in self.probe_ews:
                        self.probe_ews[k, id(ln), a] = _ew(
                            line_synthesis(ln, self.probe_rungs[a][k][index[id(ln)]]),
                            ln.fit_window, profile, ln.lam0)
                ews = np.array([self.probe_ews[k, id(ln), a] for a in offsets])
                solution = (float("nan"), float("nan"))
                if np.nanmin(ews) <= ln.ew <= np.nanmax(ews):
                    ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
                    inv = cog_invert(offsets, ews, ln.ew, ew_err=ew_err, lam0_aa=ln.lam0)
                    solution = (self.base_x + float(inv["a"]), float(inv["sigma_a"]))
                solutions[id(ln)] = solution
            self.probe_solutions.append(solutions)


def screen_lines(context, *, blend_max=0.30, weak_frac=0.10):
    """Screen every included line of ``context`` at the request's composition, in place.

    The screens use the request's solved departures when it carries them, which costs no
    solve, and LTE otherwise. They are independent of which transitions the rungs
    synthesise: a line whose element-only EW is below ``weak_frac`` of the measured EW is
    flagged ``weak_line``, and one whose EW with the element excluded exceeds ``blend_max``
    of the EW with every species is flagged ``blended``. Both are excluded.
    """
    lines = [ln for ln in context.lines if ln.include]
    if not lines:
        return
    x0 = context.base_x
    only, full, excluded = context.synthesize([
        context.variant(x0, species="only"),
        context.variant(x0, species="all"),
        context.variant(x0, species="excluded")], lte=not context.solved_departures)
    for w, ln in enumerate(context.lines):
        if not ln.include:
            continue
        fw = ln.fit_window
        ln.ew_species = float(_ew(line_synthesis(ln, only[w]), fw, "trapz", ln.lam0))
        if not np.isfinite(ln.ew_species) or ln.ew_species < weak_frac * ln.ew:
            ln.flags, ln.include = tuple(dict.fromkeys(ln.flags + ("weak_line",))), False
            continue
        blend, _, _ = blend_fraction(line_synthesis(ln, full[w]),
                                     line_synthesis(ln, excluded[w]), window=fw, x0=ln.lam0)
        ln.blend_frac = float(blend)
        if np.isfinite(blend) and blend > blend_max:
            ln.flags, ln.include = tuple(dict.fromkeys(ln.flags + ("blended",))), False


def fit_lines_cog(context, *, limit, max_rounds, profile, report=None):
    """Adaptive curve-of-growth abundances of the included lines of ``context``, in place.

    Each line's ``a`` is its [X/Fe]; the rungs are [X/Fe] offsets from the element's
    abundance in the request's composition. Returns whether every line was bracketed.
    """
    lines = [ln for ln in context.lines if ln.include]
    if not lines:
        return True
    x0 = context.base_x
    def probe_offsets():
        # the probes' solutions are bracketed like the lines' own, so every star inverts
        # through the same complete rungs
        context.solve_probes(lines, profile)
        return [a - x0 for solutions in context.probe_solutions
                for a, _sigma in solutions.values()]

    _rungs, bracketed = adaptive_cog(
        lines, context.synthesize_offsets(lines), step=context.request.config.abundance_step,
        limit=limit, max_rounds=max_rounds, profile=profile,
        cancel_token=context.cancel_token, report=report,
        also=probe_offsets if context.probes else None)
    if context.probes:
        context.solve_probes(lines, profile)
    for ln in lines:
        ln.a = x0 + ln.a if np.isfinite(ln.a) else ln.a
        ln.cog_origin += x0
        ln.a_grid = tuple(x0 + a for a in ln.a_grid)
        ln.cog_specs = {x0 + a: spec for a, spec in ln.cog_specs.items()}
    return bracketed


def fit_lines_profile(context, target, *, limit, max_rounds, profile, report=None):
    """Direct profile-fit abundances of the included lines of ``context``, in place.

    Each line's ``a`` is its [X/Fe], fitted to ``target`` on rungs at [X/Fe] offsets from
    the element's abundance in the request's composition (:func:`fit_profiles`). Returns
    whether every line was bracketed.
    """
    lines = [ln for ln in context.lines if ln.include]
    if not lines:
        return True
    _rungs, bracketed = fit_profiles(
        lines, context.synthesize_offsets(lines), target,
        step=context.request.config.abundance_step, limit=limit, max_rounds=max_rounds,
        profile=profile, cancel_token=context.cancel_token, report=report)
    x0 = context.base_x
    for ln in lines:
        ln.a = x0 + ln.a if np.isfinite(ln.a) else ln.a
        ln.cog_origin += x0
        ln.cog_specs = {x0 + a: spec for a, spec in ln.cog_specs.items()}
    return bracketed


def fit(request, *, target, line_mask=None, measured=None, method="cog",
        background_lines=True, profile="voigt", snr=None, cont_pct=0.0, limit=1.0,
        max_rounds=CORR_MAX_ROUNDS, blend_max=0.30, weak_frac=0.10, progress=None,
        cancel_token=None):
    """Measure ``target`` in the windows of ``line_mask`` and fit the element's abundance.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest`; ``target`` a
    :class:`~tumult.spectrum.Spectrum` or a spectrum file, ``line_mask`` a
    :class:`~tumult.linemask.LineMaskSet` or a mask file naming its element. The mask
    windows are the measurement windows; synthesis windows are padded internally, and the
    configured line list identifies the dominant transition of each window. ``measured``
    instead continues a result of :func:`measure_lines` (or of an earlier fit): its lines,
    as curated, are fitted without measuring the target again, and its identity and
    annotations are kept.

    Both methods measure each window's EW (``profile``, with ``snr`` and ``cont_pct``
    errors), screen weak and blended lines at the request's composition
    (:func:`screen_lines`) and work on rungs of ``config.abundance_step`` dex within
    ``+-limit`` dex of the element's abundance in that composition, each rung one synthesis
    over every window. ``method="cog"`` inverts adaptive curves of growth
    (:func:`adaptive_cog`); ``method="profile"`` fits abundance, velocity shift and
    broadening to the target on the same rungs (:func:`fit_profiles`). With
    ``background_lines=False`` the rungs synthesise the element's transitions alone; the
    screens always compare against every species. Departures requested on the request
    (``departures.request``) make every rung an NLTE solve of its own; solved departures
    are read back by the rungs and screens of either method, one solve shared by every
    rung, which is the approximation a 3D or NLTE run on fixed departures accepts.

    Returns a completed :class:`AbundanceResult` in [X/Fe] whose ``request`` and ``target``
    record the inputs; ``request.update(abundances=result)`` carries the fitted abundance on.
    """
    import uuid
    from .synthesis.request import SynthRequest

    _check_method(method)
    if not isinstance(request, SynthRequest):
        raise TypeError("abundance.fit expects a SynthRequest")
    request.check()
    spectrum, target_record = _target_spectrum(target)
    config = request.config
    report = progress or (lambda _result: None)
    if measured is not None:
        if not measured.lines:
            raise ValueError("The measured result has no lines to fit.")
        result = deepcopy(measured)
        element = result.element
    else:
        if line_mask is None:
            raise ValueError("Pass the line_mask to measure, or a measured result.")
        mset = _line_mask_set(line_mask)
        element = mset.element
        result = measure_lines(mset, spectrum,
                               linelist=str(config.resolve(config.linelist)) if config.linelist
                               else None, profile=profile, snr=snr, cont_pct=cont_pct)
    result.request, result.target, result.method = request, target_record, method
    result.cog_settings = {"background_lines": bool(background_lines), "limit": float(limit)}
    result.teff, result.logg, result.feh = request.star.teff, request.star.logg, request.star.feh
    result.vmic = float("nan") if request.star.vmic is None else request.star.vmic
    report(result)

    context = _LatticeContext(request, result.lines, element,
                              background_lines=background_lines, limit=limit,
                              cancel_token=cancel_token)
    screen_lines(context, blend_max=blend_max, weak_frac=weak_frac)
    report(result)
    if method == "cog":
        bracketed = fit_lines_cog(context, limit=limit, max_rounds=max_rounds,
                                  profile=profile, report=lambda: report(result))
    else:
        bracketed = fit_lines_profile(context, spectrum, limit=limit, max_rounds=max_rounds,
                                      profile=profile, report=lambda: report(result))
    for ln in result.lines:
        if ln.include and np.isfinite(ln.drewda):
            finalize_line(ln, spectrum, snr, cont_pct, profile=profile)
    result.relative = True
    result.snr, result.cont_pct = snr, cont_pct
    result.summary = correction_summary(result.lines)
    result.stage = "uncert"
    result.termination = "completed" if bracketed else "round_limit"
    composition = request.abundances
    if composition is not None:
        composition = composition.converted("[X/Fe]", feh=context.feh, solar=context.solar)
    result.abundances = ({"values": dict(composition), "convention": "[X/Fe]"}
                         if composition is not None else {})
    result.provenance = _provenance(request, context, background_lines)
    result.uid = result.uid or uuid.uuid4().hex[:8]
    report(result)
    return result


def _provenance(request, context, background_lines):
    from .departures import DepartureRequest
    return {"atmosphere": request.atmosphere.to_dict(), "star": dict(request.star),
            "background_lines": bool(background_lines),
            "departures": ("per_rung" if isinstance(request.departures, DepartureRequest)
                           else "none"),
            "quadrature": context.quadrature.to_dict() if context.quadrature else None}


def _reference_atmosphere(request, reference_atmosphere):
    """The 1D atmosphere or grid the LTE reference of a correction is computed on."""
    from .synthesis.inputs import AtmosphereGrid, Snapshot, coerce_atmosphere
    from .synthesis.snapshot_grid import SNAPSHOT_FORMATS
    if reference_atmosphere is not None:
        reference = coerce_atmosphere(reference_atmosphere)
    else:
        reference = request.atmosphere
        three_d = isinstance(reference, Snapshot) or (
            isinstance(reference, AtmosphereGrid)
            and str(reference.format).lower() in SNAPSHOT_FORMATS)
        if three_d:
            raise ValueError("A 3D NLTE correction needs a 1D LTE reference; pass "
                             "reference_atmosphere=AtmosphereGrid(...) or Atmosphere(...).")
    if isinstance(reference, Snapshot):
        raise ValueError("The LTE reference atmosphere must be 1D.")
    return reference


def nlte_corrections(request, *, method="cog", background_lines=False,
                     reference_atmosphere=None, profile="trapz", ew_min=5.0, limit=CORR_LIMIT,
                     max_rounds=CORR_MAX_ROUNDS, min_fit_points=MIN_FIT_POINTS,
                     progress=None, cancel_token=None):
    """Equilibrium NLTE abundance corrections ``A(NLTE) - A(LTE)`` of a model atom's lines.

    ``request`` carries the atom as requested departures (``departures.request``) and the
    NLTE atmosphere. The atom's wavelength quadrature is prepared once, at the request's
    star and composition, and shared by every synthesis of the operation. The LTE reference
    is synthesised on ``reference_atmosphere`` (the request's atmosphere when that is 1D)
    at the request's composition, and its atom lines stronger than ``ew_min`` mÅ are
    measured over the quadrature's line windows (:func:`measure_reference_lines`). The
    NLTE curves of growth then bracket each line's reference EW on the abundance lattice,
    every rung a departure solve of its own (:func:`run_nlte_corrections`).

    ``background_lines=False`` synthesises the atom's transitions alone; with ``True`` the
    configured line list is added to the reference and every rung.

    Returns an :class:`AbundanceResult` at stage ``"corr"``; ``termination`` is
    ``"completed"`` when every line was bracketed and ``"round_limit"`` otherwise.
    """
    from .departures import DepartureRequest
    from .synthesis.request import SynthRequest

    if method != "cog":
        _check_method(method)
        raise ValueError("NLTE corrections are inverted through curves of growth; use "
                         "method='cog'.")
    if not isinstance(request, SynthRequest):
        raise TypeError("nlte_corrections expects a SynthRequest")
    departures = request.departures
    if not isinstance(departures, DepartureRequest):
        raise ValueError("NLTE corrections need the model atom as requested departures: "
                         "request.update(departures=departures.request(request, atom=...)).")
    if request.windows:
        raise ValueError("NLTE corrections measure the atom's own line windows; clear the "
                         "request's windows.")
    reference = _reference_atmosphere(request, reference_atmosphere)
    request.check()
    lte_request = request.update(atmosphere=reference, departures=None)
    lte_request.check()
    report = progress or (lambda _result: None)
    opts = CorrectionOpts(profile=profile, ew_min=ew_min, limit=limit, max_rounds=max_rounds,
                          min_fit_points=min_fit_points, line_grid=departures.line_grid,
                          linelist=background_lines)
    quadrature, quadrature_atmosphere = departures.prepare(request)
    atom = request.config.resolve(departures.atom)
    result = measure_reference_lines(lte_request, atom, opts=opts, quadrature=quadrature,
                                     cancel_token=cancel_token)
    result.request, result.method = request, "cog"
    report(result)
    if cancel_token is not None:
        cancel_token.check()
    run_nlte_corrections(request, result, atom, opts=opts, progress=progress,
                         cancel_token=cancel_token)
    result.request = request
    result.provenance = {
        "atmosphere": request.atmosphere.to_dict(), "reference_atmosphere": reference.to_dict(),
        "star": dict(request.star), "background_lines": bool(background_lines),
        "departures": "per_rung", "convention": "A(NLTE) - A(LTE)",
        "quadrature": quadrature.to_dict() if quadrature is not None else None,
        "quadrature_atmosphere": quadrature_atmosphere}
    return result


def measurement_spectrum(spec, line):
    """Temporary normalised view for EW integration; never changes the target spectrum."""
    if not line.renorm:
        return spec
    from .continuum import model_continuum
    cont = model_continuum(line.continuum, spec.xdata.to_value(u.AA))
    if np.any(cont <= 0):
        raise ValueError("Local continuum levels must be positive")
    return spec / cont


def measured_profile_func(line, target, profile):
    """The profile ``f(λ_Å) -> normalised flux`` measured for ``line`` on ``target`` over
    its fit window (through its local continuum), for a line that does not carry it:
    ``fit_func`` is a live callable and does not survive serialisation."""
    target = target.convert_units(xunit=u.AA)
    return measure_ew(measurement_spectrum(target, line), window=line.fit_window,
                      profile=profile, x0=line.lam0).func


def displayed_model(line):
    """The fitted spectrum in the prepared target's fixed coordinate system."""
    if line.fit_spec is None or not line.renorm:
        return line.fit_spec
    from .continuum import model_continuum
    return line.fit_spec * model_continuum(line.continuum,
                                          line.fit_spec.xdata.to_value(u.AA))


def refresh_line(line, profile, target=None):
    """Solve an edited line again on its cached rungs, by the method it was solved with:
    :func:`refit_line_profile` (which needs the ``target``) for a profile fit, else
    :func:`refresh_line_cog`."""
    if line.method == "fit":
        if target is None:
            raise ValueError("A profile-fitted line needs its target to be fitted again")
        refit_line_profile(line, target, profile)
    else:
        refresh_line_cog(line, profile)


def _start_refresh(line):
    """Clear an edited line's refresh flags; False (flagged ``refit_required``, without an
    abundance) when it lacks the cached lattice a refresh needs. A profile overlay stays as
    a preview while the markers move."""
    line.flags = tuple(f for f in line.flags
                       if f not in ("missing_rungs", "out_of_range", "refit_required"))
    line.missing_rungs = ()
    step = line.cog_step
    if len(line.cog_specs) < 2 or not (np.isfinite(step) and step > 0):
        line.a = line.sigma_a = line.drewda = float("nan")
        line.flags += ("refit_required",)
        return False
    return True


def _missing_rungs(line, needed):
    """Record the ``needed`` rungs ``line`` has not cached; True (and no abundance) when
    there are any."""
    offsets = line.cog_specs
    line.missing_rungs = tuple(sorted({round(float(a), 10) for a in needed
                                       if not any(abs(a - b) < 1e-8 for b in offsets)}))
    if not line.missing_rungs:
        return False
    line.a = line.sigma_a = line.drewda = float("nan")
    line.fit_spec = None
    line.flags += ("missing_rungs",)
    return True


def _lattice_cell(line, a):
    """The lattice points of ``line``'s rungs that bracket ``a``."""
    k = (a - line.cog_origin) / line.cog_step
    return [line.cog_origin + np.floor(k + 1e-9) * line.cog_step,
            line.cog_origin + np.ceil(k - 1e-9) * line.cog_step]


def refit_line_profile(line, target, profile):
    """Fit a profile-fitted line again on its cached rungs, with its current local continuum.

    No synthesis runs. The rungs span the fit's whole abundance range, so a solution at
    either end is flagged ``out_of_range``. A solution in a lattice cell whose rungs are
    not cached requests them (``missing_rungs``). Either way no abundance is reported, and
    the include choice is retained.
    """
    if not _start_refresh(line):
        return
    offsets = sorted(line.cog_specs)
    specs = [line_synthesis(line, line.cog_specs[a]) for a in offsets]
    lo, hi = offsets[0], offsets[-1]
    guess = line.a if np.isfinite(line.a) and lo <= line.a <= hi else offsets[len(offsets) // 2]
    grid, best = _fit_line_profile(line, offsets, specs, target, guess=guess)
    if not lo + 1e-3 < line.a < hi - 1e-3:
        line.a = line.sigma_a = line.drewda = float("nan")
        line.ssr, line.n_fit = float("nan"), 0
        line.fit_spec = None
        line.flags += ("out_of_range",)
        return
    if _missing_rungs(line, _lattice_cell(line, line.a)):
        return
    _profile_solution(line, grid, best, offsets, specs, profile)


def refresh_line_cog(line, profile):
    """Invert an edited EW using cached spectra and flag absent nearest lattice rungs.

    Outside the cached range, request one more lattice interval. Until the nearest
    rungs exist, no abundance or uncertainty is reported and the include choice is retained.
    """
    if not _start_refresh(line):
        return
    offsets = sorted(line.cog_specs)
    ews = [_ew(line_synthesis(line, line.cog_specs[a]), line.fit_window, profile, line.lam0)
           for a in offsets]
    line.a_grid, line.ew_grid = tuple(offsets), tuple(ews)
    if not line.ew > 0:
        line.a = line.sigma_a = line.drewda = float("nan")
        line.ssr, line.n_fit = float("nan"), 0
        line.fit_spec = None
        line.flags += ("out_of_range",)
        return
    inv = cog_invert(offsets, ews, line.ew,
                     ew_err=line.ew_err if np.isfinite(line.ew_err) else 0.0,
                     lam0_aa=line.lam0)
    guess = inv["a"]
    step = line.cog_step
    if inv["out_of_range"]:
        # Expand by one lattice interval; do not trust a far extrapolation on saturation.
        needed = [offsets[0] - step if line.ew < min(ews) else offsets[-1] + step]
    else:
        needed = _lattice_cell(line, guess)
    if _missing_rungs(line, needed):
        return
    line.a, line.sigma_a = float(guess), float(inv["sigma_a"])
    slope = inv["dadrew"]
    line.drewda = 1 / slope if np.isfinite(slope) and slope != 0 else float("nan")
    line.sigma_a_phot, line.sigma_a_cont = line.sigma_a, 0.0
    line.fit_spec = _solution_spec(line)


def _local_continuum_edit(result, line, target, *, enabled=None, nodes=None):
    """Return a remeasured copy of ``line`` with the requested local continuum.

    Enabling seeds two synth-edge nodes from the upper 10% of the prepared target.
    Explicit nodes stay fixed until edited again; no optimisation runs during an edit.
    """
    from .continuum import model_continuum, seed_local
    target = target.convert_units(xunit=u.AA)
    continuum = nodes if nodes is not None else line.continuum
    on = line.renorm if enabled is None else bool(enabled)
    if on and not len(continuum):
        continuum = seed_local(target, line.window)
    if len(continuum):
        values = model_continuum(continuum, np.asarray(line.window))
        if np.any(values <= 0):
            raise ValueError("Local continuum levels must be positive")
    edited = deepcopy(line)
    edited.renorm, edited.continuum = on, tuple(tuple(p) for p in continuum)
    prof = measure_ew(measurement_spectrum(target, edited), window=line.fit_window,
                      profile=result.profile, snr=result.snr, cont_pct=result.cont_pct,
                      x0=line.lam0)
    edited.ew, edited.ew_err = prof.maa, prof.width_err_maa
    edited.fwhm = float(prof.fwhm) if prof.fitted else float("nan")
    edited.fit_func = prof.func
    if edited.cog_specs or np.isfinite(edited.a) or edited.fit_spec is not None:
        refresh_line(edited, result.profile, target)
        if np.isfinite(edited.a) and np.isfinite(edited.drewda):
            edited.sigma_a, edited.sigma_a_phot, edited.sigma_a_cont = sigma_a_from_rew(
                prof, edited.drewda)
    return edited


def _reaggregate(result):
    if hasattr(result, "reaggregate"):
        result.reaggregate()
    else:
        result.summary = correction_summary(result.lines)


def set_local_continuum(result, line, target, *, enabled=None, nodes=None):
    """Edit one measurement's fixed local continuum and immediately remeasure its EW."""
    edited = _local_continuum_edit(result, line, target, enabled=enabled, nodes=nodes)
    line.__dict__.update(edited.__dict__)
    _reaggregate(result)
    return line


def set_local_continua(result, edits, target):
    """Remeasure several continuum edits, then aggregate once.

    ``edits`` contains ``(line, kwargs)`` pairs. A failed measurement leaves that line
    unchanged; failures are returned as ``(line, exception)`` pairs.
    """
    changed, failed = [], []
    for line, kwargs in edits:
        try:
            changed.append((line, _local_continuum_edit(result, line, target, **kwargs)))
        except (ValueError, RuntimeError) as exc:
            failed.append((line, exc))
    for line, edited in changed:
        line.__dict__.update(edited.__dict__)
    if changed:
        _reaggregate(result)
    return failed


def complete_cog_rungs(lines, synthesize, *, profile, target=None,
                       max_rounds=CORR_MAX_ROUNDS, cancel_token=None, report=None):
    """Compute only absent rungs, solving each line again after each batch until bracketed.

    ``synthesize(abundances)`` returns spectra aligned with ``lines`` for each abundance.
    Lines are solved by their own method (:func:`refresh_line`); a profile fit needs the
    ``target``. Existing spectra are retained, and no target or continuum is modified.
    """
    pending = [ln for ln in lines if ln.missing_rungs]
    for _ in range(max_rounds):
        needed = sorted({a for ln in pending for a in ln.missing_rungs})
        if not needed:
            break
        if cancel_token is not None:
            cancel_token.check()
        batches = list(synthesize(needed))
        if len(batches) != len(needed):
            raise ValueError("Synthesis did not return every requested abundance rung")
        for a, spectra in zip(needed, batches):
            if len(spectra) != len(lines):
                raise ValueError("Synthesis did not return every line window")
            for ln, spec in zip(lines, spectra):
                if any(ln is p for p in pending):
                    ln.cog_specs.setdefault(a, spec.crop(window=ln.window))
        for ln in pending:
            refresh_line(ln, profile, target)
        if report is not None:
            report()
    return not any(ln.missing_rungs for ln in lines)


def complete_result_rungs(result, *, target=None, cancel_token=None, progress=None):
    """Complete an edited result's CoG at its stored atmosphere and composition.

    This updates line abundances and diagnostics only. A balance whose EWs changed needs
    an explicit stellar-parameter refit before it can claim convergence again.
    """
    if result.request is None:
        raise ValueError("The measurement has no stored synthesis request; refit it first")
    from .measurement import BalanceResult
    balance = isinstance(result, BalanceResult)
    request = result.request.update(star=result.params) if balance else result.request
    settings = result.cog_settings
    context = _LatticeContext(request, result.lines, result.element,
                              background_lines=settings.get("background_lines", True),
                              limit=settings.get("limit", 1.0), cancel_token=cancel_token)
    lines = context.lines
    shift = request.star.feh + result.solar_x if balance else 0.0

    def synthesize(abundances):
        rungs = [a - shift for a in abundances]
        return context.synthesize([context.variant(x) for x in rungs], rungs=rungs)

    prepared = None if target is None else target.convert_units(xunit=u.AA)
    complete_cog_rungs(lines, synthesize, profile=result.profile, target=prepared,
                       cancel_token=cancel_token,
                       report=(lambda: progress(result)) if progress else None)
    if target is not None:
        for ln in lines:
            if np.isfinite(ln.a) and np.isfinite(ln.drewda):
                finalize_line(ln, target, result.snr, result.cont_pct, result.profile)
    if balance:
        result.reaggregate()
    else:
        result.summary = correction_summary(result.lines)
    return result
