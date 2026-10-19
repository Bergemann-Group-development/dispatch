"""Intensity-array physics: the saved ``ie`` layout, the angle quadrature, and the
mu-integration to flux.

The physics-critical pieces of reading a run live here, apart from the run class:

* :func:`read_set_manifest` parses ``sets.txt``, the run's own description of what
  each slot on the saved set axis holds;
* :func:`base_set_indices` locates the three physics components in a manifest,
  for callers that split arrays other than ``ie``/``flux`` (``i3``, say);
* :func:`split_ie` replaces the historical bare ``ie.shape[-1] == 3`` check with an
  explicit layout dispatch -- a manifest lookup where the run wrote one, and the
  historical count-based guess as the fallback for runs that did not; it also splits
  the ``flux`` patch save (same component layout, ray axis already collapsed to 1);
* :func:`integrate_flux` is the disk-integration contract ``F = sum(w_mu * I_mu)``,
  used as a fallback in :meth:`tumult.run.Run.read` for runs saved before the
  Fortran side started integrating and saving ``flux`` directly;
* :func:`read_angle_quadrature` parses ``angle_quadrature.txt`` (per-ray direction
  cosines and weights; the stored ``wmu`` are normalised over the unit sphere, so the
  flux weights carry the ``2*pi`` factor).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .physics import vac2air


def read_segments(folder):
    """Parse ``save/out_segments.txt``: the windows a multi-window run synthesised.

    Written by ``spectrum_type%save`` when the output grid is segmented, one row per
    window::

        # segment  first  last  blue_vac_AA  red_vac_AA
                1      1    195    5151.43473    5171.44006
                2    196    389    5181.44272    5201.44805

    ``first``/``last`` are inclusive 1-based indices into the concatenated ``out_lam``
    grid. The concatenation has a gap at every segment boundary, so a reader cannot
    recover the windows from the wavelengths alone. Returns a list ordered by segment,
    each entry ``{'first': int, 'last': int, 'blue': float, 'red': float}`` with the
    bounds converted to air to match ``Run.lam``, or ``None`` for a single-window run,
    which has no such file and needs none.
    """
    path = Path(folder) / 'save' / 'out_segments.txt'
    if not path.exists():
        return None
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        cols = line.split()
        out.append({'first': int(cols[1]), 'last': int(cols[2]),
                    'blue': float(vac2air(np.array([float(cols[3])]))[0]),
                    'red': float(vac2air(np.array([float(cols[4])]))[0])})
    return out or None


def read_set_manifest(folder):
    """Parse ``sets.txt``: what each slot on the saved set axis actually contains.

    The file is a long table, one row per (set, element), written by
    ``rt_patch_t%save_set_manifest``::

        # set  kind       element   da
            1  continuum  -         -
            2  nlte       -         -
            3  lte        -         -

    A set with no abundance offsets is a single row with ``-`` in both the element
    and ``da`` columns; a set perturbing several elements gets one row per element,
    all sharing the same set number. Returns a list ordered by set index, each entry
    ``{'kind': str, 'offsets': {element: da}}``, or ``None`` if the run predates the
    manifest (in which case the caller falls back to the count-based layout guess).
    """
    path = Path(folder) / 'sets.txt'
    if not path.is_file():
        return None
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None,
                     names=['set', 'kind', 'element', 'da'], na_filter=False)
    if df.empty:
        raise ValueError(f"empty set manifest: {path}")
    sets = {}
    for row in df.itertuples(index=False):
        number = int(row.set)
        kind = str(row.kind).strip().lower()
        if number < 1:
            raise ValueError(f"set numbers must start at 1 in {path}")
        if kind not in {'continuum', 'nlte', 'lte'}:
            raise ValueError(f"unknown set kind {kind!r} in {path}")
        rec = sets.setdefault(number, {'kind': kind, 'offsets': {}})
        if rec['kind'] != kind:
            raise ValueError(f"set {number} has conflicting kinds in {path}")
        element = str(row.element).strip()
        if element and element != '-':
            if element in rec['offsets']:
                raise ValueError(f"set {number} repeats element {element!r} in {path}")
            try:
                rec['offsets'][element] = float(row.da)
            except ValueError as exc:
                raise ValueError(
                    f"set {number}, element {element!r} has invalid da={row.da!r}"
                ) from exc
        elif str(row.da).strip() != '-':
            raise ValueError(f"set {number} has da without an element in {path}")
    expected = list(range(1, max(sets) + 1))
    if sorted(sets) != expected:
        raise ValueError(f"set numbers must be contiguous from 1 in {path}")
    return [sets[key] for key in sorted(sets)]


def write_set_manifest(path, sets):
    """Write *sets* in the canonical long-table ``sets.txt`` format.

    ``sets`` has the structure returned by :func:`read_set_manifest`. A set with
    several abundance offsets is represented by several rows with the same number.
    The parent directory is created for synthesis cache manifests.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# set  kind       element  da\n"]
    for number, item in enumerate(sets, 1):
        kind = str(item['kind']).strip().lower()
        if kind not in {'continuum', 'nlte', 'lte'}:
            raise ValueError(f"unknown set kind {kind!r}")
        offsets = dict(item.get('offsets') or {})
        if not offsets:
            lines.append(f"{number:5d}  {kind:<9s}  {'-':<7s}  -\n")
            continue
        for element, da in offsets.items():
            element = str(element).strip()
            if not element or any(c.isspace() for c in element):
                raise ValueError(f"invalid element name {element!r}")
            lines.append(f"{number:5d}  {kind:<9s}  {element:<7s}  {float(da):+.12g}\n")
    path.write_text(''.join(lines))


def ladder_sets(sets):
    """Return ``(index, element, da)`` rows for the abundance-offset sets.

    ``index`` is the zero-based array index used by Python readers. Multi-element
    sets therefore produce several tuples sharing an index.
    """
    if not sets:
        return []
    return [(index, element, da)
            for index, item in enumerate(sets)
            for element, da in item.get('offsets', {}).items()]


def base_set_indices(sets, ncomp):
    """Locate the unperturbed physics sets in a manifest: ``(main, lte, cnt, mode)``.

    The three historical components are the sets carrying no abundance offset; ladder
    sets, when they exist, are everything else and are unpacked separately. An LTE run
    has no ``nlte`` set and its ``lte`` set *is* the intensity, which is what the
    ``mode == 'LTE'`` contract has always meant.
    """
    if len(sets) != ncomp:
        raise ValueError(
            f"set manifest lists {len(sets)} sets but the saved array has {ncomp} "
            f"components; the manifest and the run output disagree")

    def find(kind):
        for i, s in enumerate(sets):
            if s['kind'] == kind and not s['offsets']:
                return i
        return None

    i_cnt = find('continuum')
    if i_cnt is None:
        raise ValueError("set manifest has no unperturbed 'continuum' set")
    i_nlte = find('nlte')
    i_lte = find('lte')
    if i_nlte is not None:
        return i_nlte, i_lte if i_lte is not None else i_nlte, i_cnt, 'NLTE'
    if i_lte is not None:
        return i_lte, i_lte, i_cnt, 'LTE'
    raise ValueError("set manifest has no unperturbed 'nlte' or 'lte' set")


def split_ie(ie, legacy, sets=None):
    """Split a saved intensity array into ``(ie, ie_lte, ie_cnt, mode)``.

    The dispatch output stacks components on the last axis of the patch array
    ``ie[nlam, 1, 1, nang, ncomp]``. Where the run wrote a set manifest, ``sets``
    (from :func:`read_set_manifest`) says which component is which and the layout is
    a lookup. Without one the component count is all there is to go on, and which
    component is which depends on the output generation (detected upstream from the
    namelist: ``legacy`` = no ``read_departures`` key):

    * legacy (always 3 components): ``[intensity, lte, continuum]``; the run's RT mode
      is inferred by comparing the intensity to its LTE component.
    * modern: ``[intensity, continuum]`` for an LTE run (the LTE view *is* the
      intensity), ``[intensity, continuum, lte]`` for an NLTE run.

    Raises ``ValueError`` on any other component count, and on a manifest that
    disagrees with the saved array — previously an unrecognised shape was silently
    misclassified (or crashed later with an opaque IndexError).
    """
    ncomp = ie.shape[-1]
    if sets:
        i_main, i_lte, i_cnt, mode = base_set_indices(sets, ncomp)
        return (ie[:, 0, 0, :, i_main], ie[:, 0, 0, :, i_lte],
                ie[:, 0, 0, :, i_cnt], mode)
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
