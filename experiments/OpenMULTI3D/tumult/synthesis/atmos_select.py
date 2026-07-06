"""Atmosphere-node selection for stellar-grid synthesis.

Extracted from ``build_stellar_grid`` in the Phase 6 decomposition: given a target star
and an atmosphere-file grid, decide *which* model atmospheres to synthesise and with what
interpolation weights. vmic is decoupled from the spatial selection — it barely affects
the atmospheric structure but strongly affects line strengths, and MULTI3D applies it as
a uniform namelist override — so atmospheres are bracketed only in (teff, logg, feh) and
each selected node is synthesised at the ``vmic_step``-lattice points bracketing the
query vmic.
"""
from __future__ import annotations

import math

import numpy as np


def _vmic_levels(vmic, step):
    """Lattice points bracketing ``vmic`` on a ``step``-spaced grid, with weights.

    Returns ``[(v, w), ...]`` summing to 1. ``vmic`` is synthesised by overriding the
    atmosphere's microturbulence (a uniform Doppler term, cheap to vary on a fixed
    atmosphere), so we snap to a regular lattice — reused across nearby queries via the
    run cache — and interpolate linearly between the two enclosing lattice points.
    """
    vmic = max(0.0, float(vmic))
    lo = math.floor(vmic / step) * step
    hi = math.ceil(vmic / step) * step
    if math.isclose(lo, hi):
        return [(lo, 1.0)]
    alpha = (hi - vmic) / (hi - lo)
    return [(lo, alpha), (hi, 1.0 - alpha)]


def _vmic_window(vmic, step, n_points):
    """Sorted lattice values for a pchip vmic collapse: the two bracketing points
    extended outward to ``n_points``, preferring the side the query sits closer to
    (that is where the curvature matters) and never going below 0.

    On-lattice queries return the single exact point — no extra syntheses when none
    are needed.
    """
    vmic = max(0.0, float(vmic))
    lo = math.floor(vmic / step) * step
    hi = math.ceil(vmic / step) * step
    if math.isclose(lo, hi):
        return [lo]
    vals = [lo, hi]
    while len(vals) < n_points:
        below = vals[0] - step
        prefer_below = (vmic - vals[0]) <= (vals[-1] - vmic)
        if below >= -1e-9 and prefer_below:
            vals.insert(0, max(0.0, below))
        else:
            vals.append(vals[-1] + step)
    return vals


def select_nodes(star, atmos_files, method='Delaunay', vmic_step=0.5, absmet_files=None,
                 vmic_points=2):
    """Select the synthesis nodes enclosing ``star``: atmospheres x vmic lattice points.

    Collapses the atmosphere grid to unique (teff, logg, feh) — for each spatial node the
    representative atmosphere is the one whose grid vmic is closest to the query (its
    structure best matches; vmic is overridden anyway; atmos keys may be 3- or 4-tuples) —
    selects the enclosing spatial nodes via ``get_interp_weights`` (``method`` =
    'Delaunay' or 'NaturalNeighbor'), and crosses them with the vmic lattice from
    :func:`_vmic_levels`. A node's combined weight is the spatial node weight times the
    vmic interpolation weight; the same atmosphere may appear at two vmic levels
    (distinct keys).

    ``vmic_points`` > 2 requests that many vmic lattice points (via
    :func:`_vmic_window`) for a monotone-cubic (pchip) collapse instead of the 2-point
    linear bracket. pchip weights depend on the data, so in that mode each node's
    ``weight`` is the *spatial* weight alone — the vmic axis must be collapsed by the
    consumer (``build_stellar_grid`` marks it separable on the returned grid).

    Returns
    -------
    nodes : list of (jkey, teff, logg, feh, vmic, weight, atmos_path, absmet_file)
        with ``jkey = ((teff, logg, feh), vmic)`` the job key.
    vmic : float
        The query vmic (0.0 when ``star`` has no fourth component).
    """
    from tumult.spectrum_grid import get_interp_weights

    teff, logg, feh = star[:3]
    vmic = star[3] if len(star) == 4 else 0.0

    rep = {}  # (teff, logg, feh) -> (grid_vmic, path)
    for key, path in atmos_files.items():
        spatial = tuple(float(x) for x in key[:3])
        kvmic = float(key[3]) if len(key) == 4 else 0.0
        cur = rep.get(spatial)
        if cur is None or abs(kvmic - vmic) < abs(cur[0] - vmic):
            rep[spatial] = (kvmic, path)

    spatial_keys = list(rep.keys())
    points = np.array(spatial_keys, dtype=float)
    query = np.array([teff, logg, feh])

    # select interpolation nodes (only in the dimensions where the query is off-grid)
    weights, indices = get_interp_weights(points, query, method)

    if vmic_points > 2:
        levels = [(v, None) for v in _vmic_window(vmic, vmic_step, vmic_points)]
    else:
        levels = _vmic_levels(vmic, vmic_step)

    nodes = []
    for w_spatial, idx in zip(weights, indices):
        spatial = spatial_keys[idx]
        teff_i, logg_i, feh_i = spatial
        atmos_path = rep[spatial][1]
        absmet_file = absmet_files.get(feh_i, '') if absmet_files else ''
        for v, w_vmic in levels:
            weight = float(w_spatial) if w_vmic is None else float(w_spatial) * float(w_vmic)
            nodes.append(((spatial, v), teff_i, logg_i, feh_i, v,
                          weight, atmos_path, absmet_file))
    return nodes, vmic
