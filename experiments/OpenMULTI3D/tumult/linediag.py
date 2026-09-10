"""Diagnostic measurement of line-mask windows — the shared core behind auto-generation
and the abundance / ionisation-balance measurement step.

Both do the same thing to a spectrum, differing only in where the windows come from:

* **auto-generation** (:func:`generate_masks`) starts from candidate rows of a line list and
  *derives* each window from the reference itself (``abundance.measure_line``'s inflection
  auto-crop), producing the masks of a new :class:`linemask.LineMaskSet`;
* **measuring** (:func:`measure_masks`) takes windows that already exist — hand-drawn,
  imported from a file or just generated — and measures them.

What a run measures is an artifact of one (spectrum, line list, settings) triple, so it never
touches the masks: :func:`measure_masks` returns the :class:`abundance.LineResult` list the
wizards consume, and the mask set stays what the human decided on — windows, padding and the
element they diagnose.
"""
from __future__ import annotations

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


def _line(lam0, mask, *, ion="", e_low=float("nan"), loggf=float("nan"),
          ew=float("nan"), ew_err=float("nan"), fwhm=float("nan"), flags=()):
    """One measured line: the atomic identity, the EW and its error, the FWHM and the
    screening flags, on a copy of the window they were measured in."""
    return ab.LineResult(
        lam0=float(lam0), ion=str(ion), e_low=float(e_low), loggf=float(loggf),
        ew=float(ew), ew_err=float(ew_err), fwhm=float(fwhm),
        mask=LineMask(lam0=float(lam0), lmin=mask.lmin, lmax=mask.lmax,
                      smin=mask.smin, smax=mask.smax, note=mask.note),
        flags=tuple(flags), include=True)


def _measure(ref_spec, lam0, fit_window, *, profile, snr, cont_pct):
    """EW error terms for a fixed window; ``(ew_err, ew_err_phot)``, NaN when unmeasurable."""
    try:
        prof = ab.measure_ew(ref_spec, window=fit_window, snr=snr, cont_pct=cont_pct,
                             profile=profile, x0=lam0)
        return prof.width_err_maa, prof.err_phot_maa
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
    core, off-centre fit, EW error above the EW) are kept as *switched-off* masks, so the user
    can see and override the rejection. Returns ``list[LineMask]``: the measurement decides
    which windows survive, and is then discarded — :func:`measure_masks` is what produces
    numbers a caller keeps.

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
        common = dict(lam0=lam0, lmin=fit_window[0], lmax=fit_window[1],
                      smin=fit_window[0] - pad, smax=fit_window[1] + pad, note=str(row.ion))
        if "no_window" in flags or lprof is None:
            masks.append(LineMask(**common, include=False))   # no fittable core
            continue
        ew = lprof.maa
        if not np.isfinite(ew):
            continue
        if not flags and ew < ew_min:
            continue                 # an accepted line below the EW floor
        _ew_err, ew_err_phot = _measure(ref_spec, lam0, fit_window,
                                        profile=profile, snr=snr, cont_pct=cont_pct)
        flags = _ew_unc_flag(flags, ew, ew_err_phot)
        masks.append(LineMask(**common, include=not flags))
    return masks


def measure_masks(ref_spec, masks, *, profile, snr, cont_pct, catalog=None, element=""):
    """Measure existing windows on ``ref_spec``; ``list[abundance.LineResult]``, one per mask.

    Windows are taken as given — nothing is re-derived, added or dropped, so a curated set
    keeps exactly the masks the user chose, and the masks themselves are left untouched.
    ``catalog`` supplies the atomic identity of whichever transition dominates each window;
    without it the identity stays blank and only the EW is measured. Pass ``element``
    alongside a catalog covering *every* species (not one filtered to the element) to have
    windows another species dominates flagged — see :func:`_identify`.
    """
    lines = []
    for m in masks:
        if not (np.isfinite(m.lmin) and np.isfinite(m.lmax)):
            lam0 = m.lam0 if np.isfinite(m.lam0) else 0.0
            lines.append(_line(lam0, m, flags=("no_window",)))
            continue
        lam0 = m.lam0 if np.isfinite(m.lam0) else 0.5 * (m.lmin + m.lmax)
        ident, flags = _identify(catalog, m.lmin, m.lmax, element)
        ew = fwhm = float("nan")
        try:
            prof = ab.measure_ew(ref_spec, window=(m.lmin, m.lmax), snr=snr,
                                 cont_pct=cont_pct, profile=profile, x0=lam0)
            ew, fwhm = prof.maa, float(prof.fwhm)
            ew_err, ew_err_phot = prof.width_err_maa, prof.err_phot_maa
        except Exception:
            ew_err = ew_err_phot = float("nan")
            flags = flags + ("no_window",)
        flags = _ew_unc_flag(flags, ew, ew_err_phot)
        lines.append(_line(lam0, m, ew=ew, ew_err=ew_err, fwhm=fwhm, flags=flags, **ident))
    return lines


def _identify(catalog, lmin, lmax, element=""):
    """Atomic identity of the strongest catalog transition inside ``[lmin, lmax]``.

    A window generally spans several transitions, so this is explicitly the *dominant* one and
    lives only in the measured line — never on the mask itself. Disentangling a real
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
    """Generate a :class:`~tumult.linemask.LineMaskSet` of windows for ``element``.

    Pulls every transition of ``element`` from ``linelist`` that falls inside ``ref_spec``'s
    range and derives a window for each (:func:`generate_masks`). What comes back is a set of
    windows like any other; measuring them is :func:`measure_masks`, the wizards' step 1.

    No synthesis and no atmosphere — this needs nothing but the spectrum and the list. Raises
    :class:`ValueError` when the selection yields nothing.
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

    masks = generate_masks(ref_spec, cand, pad=pad, profile=profile, ew_min=float(ew_min),
                           snr=snr, cont_pct=cont_pct)
    if not masks:
        raise ValueError(f"No {element} lines above {float(ew_min):.0f} mÅ in {reference_name}.")

    return LineMaskSet(masks=masks, pad=pad, element=element,
                       label=f"{element} · auto ({reference_name})")
