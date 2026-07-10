"""Intensity-array physics: the saved ``ie`` layout, the angle quadrature, and the
mu-integration to flux.

Extracted from the run class in the Phase 6 decomposition so the physics-critical
pieces are named, documented and unit-tested instead of living inline:

* :func:`split_ie` replaces the historical bare ``ie.shape[-1] == 3`` check with an
  explicit layout dispatch that fails loudly on an unrecognised shape; it also splits
  the ``flux`` patch save (same component layout, ray axis already collapsed to 1);
* :func:`integrate_flux` is the disk-integration contract ``F = sum(w_mu * I_mu)``,
  used as a fallback in :meth:`tumult.run.tumult.read` for runs saved before the
  Fortran side started integrating and saving ``flux`` directly;
* :func:`read_angle_quadrature` parses ``angle_quadrature.txt`` (per-ray direction
  cosines and weights; the stored ``wmu`` are normalised over the unit sphere, so the
  flux weights carry the ``2*pi`` factor).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def split_ie(ie, legacy):
    """Split a saved intensity array into ``(ie, ie_lte, ie_cnt, mode)``.

    The dispatch output stacks components on the last axis of the patch array
    ``ie[nlam, 1, 1, nang, ncomp]``; which component is which depends on the output
    generation (detected upstream from the namelist: ``legacy`` = no
    ``read_departures`` key):

    * legacy (always 3 components): ``[intensity, lte, continuum]``; the run's RT mode
      is inferred by comparing the intensity to its LTE component.
    * modern: ``[intensity, continuum]`` for an LTE run (the LTE view *is* the
      intensity), ``[intensity, continuum, lte]`` for an NLTE run.

    Raises ``ValueError`` on any other component count — previously an unrecognised
    shape was silently misclassified (or crashed later with an opaque IndexError).
    """
    ncomp = ie.shape[-1]
    if legacy:
        if ncomp != 3:
            raise ValueError(
                f"legacy dispatch output must save 3 intensity components "
                f"[ie, lte, cnt]; got {ncomp} (shape {ie.shape})")
        inten = ie[:, 0, 0, :, 0]
        inten_lte = ie[:, 0, 0, :, 1]
        inten_cnt = ie[:, 0, 0, :, 2]
        mode = 'LTE' if np.all(inten_lte == inten) else 'NLTE'
        return inten, inten_lte, inten_cnt, mode
    if ncomp == 3:
        inten = ie[:, 0, 0, :, 0]
        inten_cnt = ie[:, 0, 0, :, 1]
        inten_lte = ie[:, 0, 0, :, 2]
        return inten, inten_lte, inten_cnt, 'NLTE'
    if ncomp == 2:
        inten = ie[:, 0, 0, :, 0]
        inten_cnt = ie[:, 0, 0, :, 1]
        return inten, inten, inten_cnt, 'LTE'
    raise ValueError(
        f"unrecognised intensity layout: expected 2 (LTE) or 3 (NLTE) components on "
        f"the last axis, got {ncomp} (shape {ie.shape})")


def integrate_flux(ie, wts):
    """Disk-integrated flux from per-ray intensities: ``F(lam) = sum_mu w_mu I_mu``.

    ``wts`` are the quadrature weights from :func:`read_angle_quadrature` (they carry
    the ``2*pi`` steradian factor); the rays live on the last axis of ``ie``.
    """
    return np.sum(wts * ie, axis=-1)


def read_angle_quadrature(file):
    """Parse ``angle_quadrature.txt``: per-ray weights and direction cosines.

    Returns a dict with ``vec`` (the (mux, muy, muz) rows), ``wts`` (the flux weights,
    ``wmu * 2*pi``), ``n_ang``, and the distinct outgoing ``mus`` (mu = -muz, sorted
    ascending by |mu| convention preserved from the run class) with ``imus`` mapping
    each ray to its mu group and ``nmu`` their count. ``None`` if the file is absent.
    """
    if not file.exists():
        return None
    df = pd.read_csv(file, header=None, skiprows=1, sep=r'\s+',
                     names=['ang', 'wmu', 'mux', 'muy', 'muz'])
    vec = df.values[:, 2:]
    wts = df.values[:, 1] * 2 * np.pi
    mus, imus = np.unique(-vec[:, -1], return_inverse=True)
    return {"vec": vec, "wts": wts, "n_ang": wts.size,
            "mus": -mus, "imus": imus, "nmu": mus.size}
