"""Diagnostic measurement of line-mask windows — the shared core behind auto-generation
and the abundance / ionisation-balance measurement step.

Both do the same thing to a spectrum, differing only in where the windows come from:

* **auto-generation** (:func:`generate_masks`) starts from candidate rows of a line list and
  *derives* each window from the reference itself (``abundance.measure_line``'s inflection
  auto-crop), producing the masks of a new :class:`linemask.LineMaskSet`;
* **re-measuring** (:func:`measure_masks`) takes windows that already exist — hand-drawn or
  imported from a file — and only measures them.

Either way each mask ends up carrying a ``diag`` row: the atomic identity, the EW and its
error, the FWHM and the screening flags. That row is an artifact of one (spectrum, line list,
settings) run, so it lives on the mask in memory and in session saves but is never written to
the line-mask file — the file holds the windows a human decided on, not what one particular
observation said about them. :func:`line_results` turns a measured set back into the
:class:`abundance.LineResult` list the wizards consume.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import abundance as ab
from .linelist import strip_ion_stages
from .linemask import DEFAULT_PAD, LineMask

# the reason a candidate was screened out; a flagged mask is created but switched off, so the
# user sees why rather than wondering where the line went
FLAG_LABELS = {
    "no_window": "no fittable core",
    "off_center": "fit centre falls outside the window",
    "ew_unc": "EW error exceeds the EW",
    "no_lines": "no line-list match in the window",
    "off_element": "dominated by another species",
}


def _diag(*, ion="", e_low=float("nan"), loggf=float("nan"), ew=float("nan"),
          ew_err=float("nan"), fwhm=float("nan"), flags=()):
    """One measurement row, as stored on ``LineMask.diag`` (plain data — it round-trips
    through the session snapshot as JSON)."""
    return {"ion": str(ion), "e_low": float(e_low), "loggf": float(loggf),
            "ew": float(ew), "ew_err": float(ew_err), "fwhm": float(fwhm),
            "flags": list(flags)}


def _measure(ref_spec, lam0, fit_window, *, profile, snr, cont_pct):
    """EW error terms for a fixed window; ``(ew_err, ew_err_phot)``, NaN when unmeasurable."""
    try:
        prof = ab.measure_ew(ref_spec, window=fit_window, snr=snr, cont_pct=cont_pct,
                             profile=profile, x0=lam0)
        return prof.width_err_mAA, prof.err_phot_mAA
    except Exception:
        return float("nan"), float("nan")


def _ew_unc_flag(flags, ew, ew_err_phot):
    """Flag a line whose photon EW error exceeds its own EW (fractional uncertainty > 1)."""
    if np.isfinite(ew_err_phot) and ew > 0 and ew_err_phot / ew > 1.0:
        return tuple(flags) + ("ew_unc",)
    return tuple(flags)


def generate_masks(ref_spec, cand, *, profile, ew_min, snr, cont_pct, pad=DEFAULT_PAD):
    """Derive a window per candidate line-list row and measure it on ``ref_spec``.

    ``cand`` needs ``lam0`` / ``ion`` / ``evi`` / ``loggf`` columns. Hard fit failures and
    clean lines below ``ew_min`` are dropped; lines screened by a sanity check (no fittable
    core, off-centre fit, EW error above the EW) are kept as *switched-off* masks carrying
    the flag, so the user can see and override the rejection. Returns ``list[LineMask]``.

    Only the *fit* window is derived from the spectrum; the synth window is the fit window
    widened by ``pad`` Å, exactly as for a drawn or imported mask — so the sidebar's synth
    padding governs every set the same way and keeps governing it after generation.
    """
    masks = []
    for row in cand.itertuples():
        lam0 = float(row.lam0)
        try:
            lprof, _window, fit_window, flags = ab.measure_line(
                ref_spec, {"lam0": lam0}, profile=profile)
        except Exception:
            continue                 # hard fit failure (too few points) — drop it
        ion = str(row.ion)
        common = dict(lam0=lam0, lmin=fit_window[0], lmax=fit_window[1],
                      smin=fit_window[0] - pad, smax=fit_window[1] + pad, note=ion)
        if "no_window" in flags or lprof is None:
            # no fittable core: keep the window, flagged and off, with no EW measured
            masks.append(LineMask(**common, include=False, diag=_diag(
                ion=ion, e_low=float(row.evi), loggf=float(row.loggf), flags=flags)))
            continue
        ew, fwhm = lprof.mAA, float(lprof.fwhm)
        if not np.isfinite(ew):
            continue
        if not flags and ew < ew_min:
            continue                 # an accepted line below the EW floor
        ew_err, ew_err_phot = _measure(ref_spec, lam0, fit_window,
                                       profile=profile, snr=snr, cont_pct=cont_pct)
        flags = _ew_unc_flag(flags, ew, ew_err_phot)
        masks.append(LineMask(**common, include=not flags, diag=_diag(
            ion=ion, e_low=float(row.evi), loggf=float(row.loggf),
            ew=ew, ew_err=ew_err, fwhm=fwhm, flags=flags)))
    return masks


def measure_masks(ref_spec, masks, *, profile, snr, cont_pct, catalog=None, element=""):
    """Measure *existing* windows on ``ref_spec``, filling each mask's ``diag`` in place.

    Windows are taken as given — nothing is re-derived, added or dropped, so a curated set
    keeps exactly the masks the user chose. ``catalog`` supplies the atomic identity of
    whichever transition dominates each window; without it the identity stays blank and only
    the EW is measured. Pass ``element`` alongside a catalog covering *every* species (not
    one filtered to the element) to have windows another species dominates flagged — see
    :func:`_identify`.
    """
    for m in masks:
        fit_window = (m.lmin, m.lmax)
        if not (np.isfinite(m.lmin) and np.isfinite(m.lmax)):
            m.diag = _diag(flags=("no_window",))
            continue
        lam0 = m.lam0 if np.isfinite(m.lam0) else 0.5 * (m.lmin + m.lmax)
        ident, flags = _identify(catalog, m.lmin, m.lmax, element)
        ew = fwhm = float("nan")
        try:
            prof = ab.measure_ew(ref_spec, window=fit_window, snr=snr, cont_pct=cont_pct,
                                 profile=profile, x0=lam0)
            ew, fwhm = prof.mAA, float(prof.fwhm)
            ew_err, ew_err_phot = prof.width_err_mAA, prof.err_phot_mAA
        except Exception:
            ew_err = ew_err_phot = float("nan")
            flags = flags + ("no_window",)
        flags = _ew_unc_flag(flags, ew, ew_err_phot)
        m.diag = _diag(ew=ew, ew_err=ew_err, fwhm=fwhm, flags=flags, **ident)
    return masks


def _identify(catalog, lmin, lmax, element=""):
    """Atomic identity of the strongest catalog transition inside ``[lmin, lmax]``.

    A window generally spans several transitions, so this is explicitly the *dominant* one and
    lives only in the measurement artifact — never on the mask itself. Disentangling a real
    blend is the job of a synthesis run's ``linelist_out.txt``, which weights every transition
    by its actual contribution.

    Which transition dominates depends on the star, so what ``catalog`` carries decides how
    much this is worth: a frame with a ``strength`` column computed on a model atmosphere
    (:func:`linestrength.compute_band`) ranks by the line's actual line-to-continuum
    extinction there, while a bare line list can only fall back on the first row in the
    window. ``element``, when given, is the species the windows are *supposed* to diagnose:
    a window whose dominant transition belongs to something else is flagged rather than
    silently attributed, since the caller (ionization balance, an abundance) is about to
    read a stage or an abundance off that identity.
    """
    if catalog is None or not len(catalog):
        return {}, ()
    sub = catalog[(catalog.lam0 >= lmin) & (catalog.lam0 <= lmax)]
    if not len(sub):
        return {}, ("no_lines",)
    row = sub.loc[sub.strength.idxmax()] if "strength" in sub.columns else sub.iloc[0]
    ion = str(row.ion)
    flags = () if (not element or ion.split()[:1] == [str(element)]) else ("off_element",)
    return ({"ion": ion, "e_low": float(row.evi), "loggf": float(row.loggf)}, flags)


def line_from_mask(m, fallback_ion=""):
    """One measured mask as an :class:`abundance.LineResult`.

    Copies the window and everything the diagnostic run wrote onto the mask, so the result
    stops depending on the mask object: the set stays free to be edited afterwards without
    changing what was measured. ``fallback_ion`` names the species when the run could not
    identify one (the set's element, where the caller is content with that).
    """
    d = m.diag or {}
    lam0 = m.lam0 if np.isfinite(m.lam0) else 0.5 * (m.lmin + m.lmax)
    return ab.LineResult(
        lam0=float(lam0), ion=str(d.get("ion") or fallback_ion),
        e_low=float(d.get("e_low", float("nan"))),
        loggf=float(d.get("loggf", float("nan"))),
        ew=float(d.get("ew", float("nan"))),
        ew_err=float(d.get("ew_err", float("nan"))),
        fwhm=float(d.get("fwhm", float("nan"))),
        mask=LineMask(lam0=float(lam0), lmin=m.lmin, lmax=m.lmax,
                      smin=m.smin, smax=m.smax, note=m.note),
        flags=tuple(d.get("flags") or ()), include=True)


def line_results(mset, *, included_only=True):
    """The measured set as ``list[abundance.LineResult]``.

    Each kept mask becomes one line, carrying its windows and whatever the diagnostic run
    measured. Masks the user switched off are dropped by default: the toggle is the curation.
    """
    return [line_from_mask(m, mset.element) for m in mset.masks
            if m.include or not included_only]


def spectrum_range_aa(spec):
    """The spectrum's wavelength span as ``(lo, hi)`` in Angstrom, ascending.

    A unit flip (wavelength -> energy) reverses the order, so the result is sorted rather
    than assumed: candidate selection bisects on it.
    """
    import astropy.units as u

    x, _y = spec.values
    lo = (float(x.min()) * spec.xunit).to_value(u.AA, equivalencies=u.spectral())
    hi = (float(x.max()) * spec.xunit).to_value(u.AA, equivalencies=u.spectral())
    return (lo, hi) if lo <= hi else (hi, lo)


def select_candidates(master, element, lo, hi):
    """Rows of ``master`` for ``element`` with ``lam0`` inside ``[lo, hi]``, sorted by lam0.

    The line list names a *species* ("Fe I", "Fe II", "CH"); a diagnostic set is per element,
    so the ionisation stage is dropped and every stage of the element is a candidate.
    Molecular species never match an element symbol and so are excluded.
    """
    cand = master[(strip_ion_stages(master.ion) == element)
                  & (master.lam0 >= lo) & (master.lam0 <= hi)]
    return cand.sort_values("lam0")


def generate_mask_set(ref_spec, element, *, linelist, reference_name="spectrum",
                      profile="gauss", ew_min=5.0, snr=None, cont_pct=0.0,
                      pad=DEFAULT_PAD):
    """Generate a measured :class:`~tumult.linemask.LineMaskSet` for ``element``.

    Pulls every transition of ``element`` from ``linelist`` that falls inside ``ref_spec``'s
    range, derives and measures a window for each (:func:`generate_masks`), and records what
    produced the set in ``diag_meta`` so a later consumer can tell whether the diagnostics
    still describe the spectrum in front of it.

    No synthesis and no atmosphere — this is the wizards' step 1, needing nothing but the
    spectrum and the list. Raises :class:`ValueError` when the selection yields nothing.
    """
    from .linelist import load_master
    from .linemask import LineMaskSet, is_element

    element = str(element or "").strip().capitalize()
    if not is_element(element):
        raise ValueError(f"“{element}” is not an element symbol.")
    lo, hi = spectrum_range_aa(ref_spec)
    cand = select_candidates(load_master(linelist), element, lo, hi)
    if cand.empty:
        raise ValueError(f"No {element} lines in {lo:.1f}–{hi:.1f} Å of {reference_name}.")

    settings = dict(profile=profile, ew_min=float(ew_min), snr=snr, cont_pct=cont_pct)
    masks = generate_masks(ref_spec, cand, pad=pad, **settings)
    if not masks:
        raise ValueError(f"No {element} lines above {float(ew_min):.0f} mÅ in {reference_name}.")

    mset = LineMaskSet(masks=masks, pad=pad, element=element,
                       label=f"{element} · auto ({reference_name})")
    mset.diag_meta = {"reference": reference_name, "element": element,
                      "linelist": Path(linelist).name,
                      **{k: v for k, v in settings.items() if v is not None}}
    return mset
