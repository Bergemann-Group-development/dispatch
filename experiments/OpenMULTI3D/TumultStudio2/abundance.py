"""Mask-free, uncertainty-aware abundance pipeline (TumultStudio2 layer).

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

from dataclasses import dataclass, fields

import astropy.units as u
import numpy as np

from tumult import Spectrum
from tumult.abundance import (  # noqa: F401 (re-exports; numerics live in the library)
    aggregate,
    cog_invert,
    profile_fwhm,
    reduced_ew,
    sigma_a_from_rew,
    velocity_hwhm,
)
from tumult.profile import get_profile_class
from tumult.utils import auto_crop_core

# the per-line synthesis / jump window spans SYNTH_HWHM_MULT half-widths (HWHM) either side
# of the *fitted* line centre, but is always widened to reach at least FIT_MARGIN_HWHM HWHM
# beyond the auto-crop fit window, so the fit region stays fully enclosed. This window sets
# both the jump view range and what gets synthesised, so the preview matches the synthesis.
SYNTH_HWHM_MULT = 3.0
FIT_MARGIN_HWHM = 1.0


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
    window: tuple = (float("nan"), float("nan"))  # (xmin, xmax) [AA] — jump & synth range
    fit_window: tuple = (float("nan"), float("nan"))  # (xmin, xmax) [AA] — auto-crop fit window
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

    # --- session persistence (see AbundanceResult / the session snapshot) -------------
    _TUPLE_FIELDS = ("window", "fit_window", "a_grid", "ew_grid", "flags")

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
    """
    ew_all = _ew(spec_all, window, profile, x0)
    ew_excl = _ew(spec_excl, window, profile, x0)
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
    element: str                       # species, e.g. "Fe I"
    symbol: str                        # bare element symbol, e.g. "Fe"
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
