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
    return " ".join(t for t in (str(base).strip(), f"{element}={val:.4f}") if t)


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
    ew: float = float("nan")    # measured equivalent width [mAA]
    ew_err: float = float("nan")  # combined 1σ EW error [mAA]
    ew_species: float = float("nan")  # query-species-only EW in the fit window [mAA] (weak screen)
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
    fit_spec: object = None     # best-fit synthetic Spectrum, for the jump overlay
    fit_func: object = None     # fitted profile f(λ_Å)->flux, for the EW-area hatch on jump
    method: str = "fit"         # "fit" (direct) or "cog"
    flags: tuple = ()
    include: bool = True        # user toggle: count this line in the element mean

    # --- window views: the fit/synth windows now live on the embedded LineMask, but the
    # pipeline (synth.fit_line, measure_ew, canvas.set_abund_windows, balance/corr tools) still
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
    error terms. The returned Profile exposes the width (``.mAA``), the reduced EW (``.rew``)
    and the error split off one interface.
    """
    return spec.measure(window or None, profile=profile, x0=x0, auto_crop=False, **kw)


def _ew(spec, window=None, profile="trapz", x0=None):
    """Noiseless *width* EW of a (synthetic) spectrum over ``window``, in mÅ."""
    return measure_profile(spec, window, profile=profile, x0=x0).mAA


def _rew(spec, window=None, profile="trapz", x0=None):
    """Noiseless reduced EW ``log10(W/λ)`` of a (synthetic) spectrum over ``window`` (dimensionless)."""
    return measure_profile(spec, window, profile=profile, x0=x0).rew


def measure_ew(spec, window=None, snr=None, cont_pct=0.0, profile="trapz", x0=None):
    """Measure a line's EW as a :class:`~tumult.profile.Profile` with a photon / continuum error split.

    ``snr`` None (a synthetic reference) drops the photon term; ``cont_pct`` 0 drops the continuum
    term. A single measurement carries both shares — ``.err_phot_mAA`` (SNR) and ``.err_cont_mAA``
    (continuum offset), combining in quadrature to ``.width_err_mAA`` — so no second pass is needed.
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


def build_cog(synth_line_fn, a_ladder, window=None, profile="trapz", x0=None,
              batch_synth_fn=None):
    """Synthesise the line at each abundance in ``a_ladder`` and measure its *width* EW.

    ``synth_line_fn(a) -> Spectrum`` returns the continuum-normalised synthetic line at
    abundance ``a`` (caller's convention — [X/Fe] or A(X)). Returns ``(a_grid, ew_grid)`` (mÅ).

    ``batch_synth_fn(a_list) -> [Spectrum]`` (optional): when given, the whole ladder is
    synthesised in one call so the front-end can fan every point's run across a single process
    pool; the returned spectra must align with ``a_list``. Equivalent to the serial
    ``synth_line_fn`` loop (same per-point spectra), just scheduled concurrently — used when
    omitted.
    """
    a = np.asarray(a_ladder, dtype=float)
    if batch_synth_fn is not None:
        specs = batch_synth_fn(list(a))
        ews = np.array([_ew(s, window, profile, x0) for s in specs])
    else:
        ews = np.array([_ew(synth_line_fn(av), window, profile, x0) for av in a])
    return a, ews


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
    x, y = spec.values
    try:
        xc, _ = auto_crop_core(x, y, lam0, exclude_wing=exclude_wing, resolution=resolution)
        if len(xc) >= 3:
            return (float(np.min(xc)), float(np.max(xc))), True
    except Exception:
        pass
    return (float(lam0 - fallback), float(lam0 + fallback)), False


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
class AbundanceResult:
    """One completed element run: the per-line measurements plus the aggregated mean.

    A plain data container (no GUI deps) so it can be serialised into a session snapshot and
    re-consumed by a later Fe-balance / stellar-parameter solver. ``feh``/``relative`` record
    the abundance convention so ``a`` values are interpretable.
    """
    # Both hold the element, e.g. "Fe" — never a species like "Fe I": the abundance is per
    # element (MULTI3D scales by element, and the balance appends the stage itself). They are
    # set alike at every construction site; ``symbol`` is kept because saved sessions carry it.
    element: str
    symbol: str
    lines: list                        # list[LineResult]
    summary: dict                      # aggregate(...) output
    stage: str = "lines"               # "lines" → "synth" → "uncert" (how far it has run)
    profile: str = "gauss"             # EW/FWHM profile chosen at line-selection
    teff: float = float("nan")
    logg: float = float("nan")
    feh: float = float("nan")
    vmic: float = float("nan")
    snr: float | None = None
    cont_pct: float = float("nan")
    relative: bool = True              # True → ``a`` is [X/Fe]; False → A(X)
    uid: str = ""
    name: str = ""

    # --- session persistence ------------------------------------------------------------
    def to_dict(self):
        """A JSON-safe view for the session snapshot. Captures the per-line records plus the
        reference the tools attach dynamically: ``ref_uid`` (AbundTool → a canvas observation)
        and ``ref_spec`` (CorrTool → the synthetic 1D LTE reference, embedded)."""
        d = {}
        for f in fields(self):
            d[f.name] = [ln.to_dict() for ln in self.lines] if f.name == "lines" \
                else getattr(self, f.name)
        if hasattr(self, "ref_uid"):
            d["ref_uid"] = self.ref_uid
        ref_spec = getattr(self, "ref_spec", None)
        if ref_spec is not None:
            d["ref_spec"] = ref_spec.to_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        d = dict(d or {})
        names = {f.name for f in fields(cls)}
        lines = [LineResult.from_dict(x) for x in (d.get("lines") or [])]
        kw = {k: v for k, v in d.items() if k in names and k != "lines"}
        res = cls(lines=lines, **kw)
        if "ref_uid" in d:
            res.ref_uid = d["ref_uid"]
        if d.get("ref_spec") is not None:
            res.ref_spec = Spectrum.from_dict(d["ref_spec"])
        return res


@dataclass
class Thresholds:
    ew_min: float = 5.0         # mÅ — below this a line is too weak to trust
    logrw_max: float = -4.8     # reduced-EW ceiling — above this the line is saturating
    blend_max: float = 0.30     # contaminant EW fraction ceiling
    weak_frac: float = 0.10     # species-only EW must be ≥ this fraction of the measured EW


# --- staged pipeline (the 3-step wizard drives these one stage at a time) -----------------

def measure_line(ref_spec, line, profile="gauss", fallback=0.4, resolution=None):
    """Stage 1: fit a line on the reference and return the bound Profile plus its windows.

    Independent of any stellar parameters. Returns ``(prof, window, fit_window, flags)``:
    ``prof`` is the fitted :class:`tumult.profile.Profile` — it carries the EW (``.mAA``),
    the FWHM, the model ``.func`` (for hatching the integrated area on jump) **and the
    fit-covariance EW error** (``.width_err_mAA``; with no per-pixel noise given, the
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
    fit_window, found = auto_window(ref_spec, lam0, fallback=fallback, resolution=resolution)
    if not found:
        # no fittable core resolved: there is no trustworthy window, so don't measure an EW
        # at all — return it flagged (no Profile / NaN windows), to keep but exclude.
        nan = float("nan")
        return None, (nan, nan), (nan, nan), ("no_window",)
    sub = _crop(ref_spec, fit_window)
    x, y = sub.values
    try:
        prof = get_profile_class(profile).measure(
            np.asarray(x, float), np.asarray(y, float), x0=lam0, xunit=sub.xunit)
    except ValueError as exc:
        raise RuntimeError(f"profile fit failed at {lam0:.3f}") from exc
    center = float(prof.params[0])
    hwhm = 0.5 * float(prof.fwhm)
    # synth window: SYNTH_HWHM_MULT HWHM about the fitted centre, but at least FIT_MARGIN_HWHM
    # HWHM beyond the fit window on each side so the fit region is always enclosed.
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
             grad_step=0.125, grad_batch_fn=None):
    """Stage 2 (direct fit): screen the line, fit its abundance + vshift/vdop, and get ``dREW/dA``.

    Injected callables (synthesis-agnostic so this stays engine-independent and testable):

    * ``synth_only_fn(a) -> Spectrum`` — the line at abundance ``a`` with *only* the target
      species (``only_elements``); its fit-window width EW is the weak-line screen.
    * ``synth_all_fn(a) -> Spectrum`` — the *full* synthesis at abundance ``a`` (all species), used
      for the blend metric and the ``dREW/dA`` gradient.
    * ``synth_excl_fn() -> Spectrum`` — the window with the target species removed
      (``exclude_elements``); with ``synth_all_fn(a_guess)`` gives the blend fraction.
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
    blend, _, _ = blend_fraction(synth_all_fn(a_guess), synth_excl_fn(), window=fw,
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
    line.ew, line.ew_err = prof.mAA, prof.width_err_mAA
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

def measure_lines(mset, ref_spec, *, reference_name, linelist, profile, snr, cont_pct):
    """Stage 1: measure a mask set's switched-on windows on the reference observation.

    The windows are taken exactly as curated — nothing is re-derived, added or dropped. A set
    whose diagnostics already came from *this* reference is reused as it stands; otherwise
    the windows are re-measured here (and the set's diagnostics table follows), with the
    atomic identity of each window's dominant transition attached from the line list.

    Returns ``(result, remeasured)``; ``remeasured`` says whether ``mset`` was written to, so
    a caller with a view of the mask set knows to repaint it. Raises :class:`ValueError` when
    the set holds nothing measurable.
    """
    import uuid

    from . import linediag
    from .linelist import load_master, strip_ion_stages

    stale = (mset.diag_meta or {}).get("reference") != reference_name
    remeasured = bool(stale or not mset.has_diag)
    if remeasured:
        try:
            master = load_master(linelist)
            catalog = master[strip_ion_stages(master.ion) == mset.element]
        except Exception:
            # the catalog only names the dominant transition of each window in the
            # diagnostics artifact; it never touches the mask or the measurement
            logging.exception("Could not read the line list for %s", mset.element)
            catalog = None
        linediag.measure_masks(ref_spec, [m for m in mset.masks if m.include],
                               profile=profile, snr=snr, cont_pct=cont_pct, catalog=catalog)
        mset.diag_meta = {"reference": reference_name, "element": mset.element,
                          "linelist": Path(linelist).name, "profile": profile,
                          "snr": snr, "cont_pct": cont_pct}
    lines = linediag.line_results(mset)
    if not lines:
        raise ValueError(f"No measurable windows in “{mset.name}”.")
    result = AbundanceResult(element=mset.element, symbol=mset.element, lines=lines,
                             summary=aggregate([]), stage="lines", profile=profile,
                             uid=uuid.uuid4().hex[:8])
    result.name = f"{mset.element} · {len(lines)} lines"
    return result, remeasured


def run_abundance(request, result, observation, *, progress=None, cancel_token=None):
    """Stage 2: screen and fit every included line of ``result``, in place.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest`; ``observation`` the spectrum the
    lines were measured on. Each line gets the four syntheses :func:`fit_line` needs — the
    species in isolation, the full window, the window with the element removed, and the
    adaptive abundance fit — plus the dREW/dA flanks, which go in as one batch so they fold
    onto a single run's abundance-set axis. ``progress(result)`` is called as each line
    lands, so a caller can stream partial results; it must not raise.
    """
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

    def synth_window(rv_over, window):
        rv = dict(base_rv)
        rv.update(rv_over)
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
            return synth.compute_batch(target, prefs, rvs, *padded_window(_w), grid=grid,
                                       atmos_file=atmos_file, cancel_token=cancel_token)

        def synth_excl_fn(_w=ln.window):
            return synth_window({"exclude_elements": symbol}, _w)

        def abund_fit_fn(guess, _w=ln.window, _fw=ln.fit_window):
            return synth.fit_line(observation, target, prefs, base_rv, symbol, guess,
                                  *padded_window(_w), _fw,
                                  grid=grid, atmos_file=atmos_file,
                                  cancel_token=cancel_token)

        try:
            fit_line(ln, synth_only_fn=synth_only_fn, synth_all_fn=synth_all_fn,
                     synth_excl_fn=synth_excl_fn, abund_fit_fn=abund_fit_fn,
                     ew_meas=ln.ew, a_guess=0.0, profile=result.profile,
                     grad_batch_fn=synth_all_batch)
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
    return result


# --- NLTE abundance corrections against a 1D LTE reference ---------------------------------
# The counterpart to measure_lines / run_abundance for the case where the "observation" is a
# synthetic 1D LTE reference computed straight from a model atom. Stage 1 measures the atom's
# own lines on that reference; stage 2 builds an NLTE curve of growth per line and inverts it
# to the abundance reproducing the reference EW. The reference sits at [X/Fe]=0, so that
# abundance *is* the correction dA = A(NLTE) - A(LTE).

CORR_SPAN = 0.6        # curve-of-growth spans +/-CORR_SPAN dex in [X/Fe] about the reference
CORR_N = 7             # number of ladder points (kept small - each is an NLTE solve)
MIN_FIT_POINTS = 7     # a fit region sampled by fewer reference points than this is too sparse
                       # to trust - flag "few_points" and exclude the line by default


@dataclass
class CorrectionOpts:
    """Thresholds and ladder shape for the NLTE-correction workflow."""
    profile: str = "gauss"
    ew_min: float = 5.0                      # mAA - reference lines weaker than this are dropped
    span: float = CORR_SPAN
    points: int = CORR_N
    min_fit_points: int = MIN_FIT_POINTS


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

    Returns an :class:`AbundanceResult` at stage ``"lines"``, carrying the reference spectrum
    as ``ref_spec``. Raises :class:`ValueError` when the reference yields nothing measurable.
    """
    import uuid

    from .synthesis import workflows as synth

    opts = opts or CorrectionOpts()
    profile, ew_min = str(opts.profile), float(opts.ew_min)
    prefs, target = request.prefs, request.target

    ref_spec = synth.compute_reference(target, prefs, atom, ref_grid,
                                       cancel_token=cancel_token)
    llist = synth.read_reference_linelist(prefs)
    if llist is None or len(llist) == 0:
        raise ValueError("No linelist_out.txt was produced by the reference run.")
    cand = llist.sort_values("lam0")
    if cand.empty:
        raise ValueError("The reference run produced no atom lines.")

    lines = []
    for row in cand.itertuples():
        species, _symbol = species_symbol(getattr(row, "ion", ""))
        try:
            lprof, win, fit_window, flags = measure_line(
                ref_spec, {"lam0": float(row.lam0)}, profile=profile)
        except Exception:
            continue
        if "no_window" in flags or lprof is None or not np.isfinite(lprof.mAA):
            continue
        ew, fwhm, func = lprof.mAA, float(lprof.fwhm), lprof.func
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
        ew_err = float(lprof.width_err_mAA)
        # flag (and exclude) a line whose EW is more uncertain than the EW itself: the
        # fit-covariance error is the only EW uncertainty on the synthetic reference, so
        # ew_err > ew means the measurement is not trustworthy.
        if np.isfinite(ew_err) and ew > 0 and ew_err > ew:
            flags.append("ew_unc")
        if not flags and ew < ew_min:
            continue
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
    result.atom_name = Path(atom).name
    result.name = result.atom_name
    return result


def run_nlte_corrections(request, result, atom, *, opts=None, progress=None,
                         cancel_token=None):
    """Stage 2: invert an NLTE curve of growth to each line's correction, in place.

    ``request`` is a :class:`~tumult.synthesis.SynthRequest` naming the *NLTE target*
    atmosphere; ``result`` the stage-1 :class:`AbundanceResult` whose lines carry the
    reference EWs. Every included line gets an NLTE ladder on the picked atmosphere (model
    atom only), which is inverted at the reference EW to give ``ln.a`` = dA.

    ``progress(result)`` is called as each line lands, so a caller can stream partial results;
    it must not raise. Returns ``result`` at stage ``"corr"`` with its summary aggregated.
    """
    from .synthesis import workflows as synth
    from .synthesis.cancel import CancelledError

    opts = opts or CorrectionOpts()
    report = progress or (lambda _result: None)
    prefs, grid = request.prefs, request.grid
    target, atmos_file = request.target, request.atmos_file
    symbol = result.symbol
    span = float(opts.span)
    ladder = np.linspace(-span, span, int(opts.points))

    def _rv(a):
        return {"mode": "NLTE", "only_atom": True, "atoms": [atom],
                "abundances": f"{symbol}={a:.4f}", "relative_abundances": True}

    # set_spectrum=False leaves spectrum_params without R/daa, so only the aa_blue/aa_red
    # window is passed: the Fortran samples the line on the model atom's own frequency grid
    # (matching the 1D LTE reference) and crops it to the synth window, rather than
    # resampling onto a uniform R grid. This is also why the ladder cannot ride a segmented
    # window_file run (compute_window_batch), which fixes a constant-R grid per window.
    def nlte_synth(a, window):
        return synth.compute(target, prefs, _rv(a), *padded_window(window), grid=grid,
                             atmos_file=atmos_file, set_spectrum=False,
                             cancel_token=cancel_token)

    # the whole abundance ladder for one line synthesised in a single call, so every point's
    # dispatch run is fanned across one shared process pool instead of being issued one after
    # another (results are identical to the serial nlte_synth loop).
    def batch_nlte_synth(a_list, window):
        return synth.compute_batch(target, prefs, [_rv(a) for a in a_list],
                                   *padded_window(window), grid=grid, atmos_file=atmos_file,
                                   set_spectrum=False, cancel_token=cancel_token)

    for ln in result.lines:
        if cancel_token is not None and cancel_token.cancelled:
            raise CancelledError()
        if not ln.include:
            continue
        try:
            a_grid, ew_grid = build_cog(
                None, ladder, window=ln.fit_window, profile=result.profile, x0=ln.lam0,
                batch_synth_fn=lambda a_list, _w=ln.window: batch_nlte_synth(a_list, _w))
            # propagate the reference EW's fit-covariance error through the local CoG slope:
            # sigma(dA) = sigma_EW * |dA/dEW|. All of it is fit/model error (no continuum
            # term on a synthetic reference).
            ew_err = float(ln.ew_err) if np.isfinite(ln.ew_err) else 0.0
            inv = cog_invert(a_grid, ew_grid, ln.ew, ew_err=ew_err)
            ln.a = float(inv["a"])
            ln.sigma_a = float(inv["sigma_a"])
            ln.sigma_a_phot, ln.sigma_a_cont = ln.sigma_a, 0.0
            ln.a_grid, ln.ew_grid = tuple(a_grid), tuple(ew_grid)
            if inv["out_of_range"]:
                ln.flags = ("out_of_range",)
                ln.include = False
            else:
                # overlay the best-fit NLTE line - synthesised at the inverted abundance - so
                # jumping to a row shows its fit against the reference, mirroring the
                # abundance tool. A failed overlay must not discard the (successful)
                # correction, so guard it separately.
                try:
                    ln.fit_spec = nlte_synth(ln.a, ln.window)
                except CancelledError:
                    raise
                except Exception:
                    logging.exception("Correction: overlay synth failed %.3f", ln.lam0)
        except CancelledError:
            raise
        except Exception:
            logging.exception("Correction: line %.3f failed", ln.lam0)
            ln.flags = ("synth_failed",)
        report(result)

    result.stage = "corr"
    result.summary = correction_summary(result.lines)
    return result
