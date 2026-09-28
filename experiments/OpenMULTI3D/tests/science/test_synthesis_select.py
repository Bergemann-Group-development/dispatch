"""Atmosphere and vmic node selection (tumult/synthesis/atmos_select.py).

The lattice bracketing, the exactness on a lattice point, the weight sum and the clamp
of a below-range vmic moved to `tests/science/test_interpolation.py` as requirement S12
of `reports/test-requirements.md`. The pchip vmic window, the representative-atmosphere
collapse and the absmet lookup remain here, and have no approved requirement of their
own yet."""
import numpy as np
import pytest

from tumult.synthesis.atmos_select import _vmic_window, select_nodes


def test_vmic_window_on_lattice_is_single_point():
    assert _vmic_window(1.0, 0.5, 4) == [1.0]


def test_vmic_window_extends_toward_query_side():
    # query 1.4 sits nearer 1.5 -> third point above, fourth below
    assert _vmic_window(1.4, 0.5, 3) == [1.0, 1.5, 2.0]
    assert _vmic_window(1.4, 0.5, 4) == [0.5, 1.0, 1.5, 2.0]
    # query 1.1 sits nearer 1.0 -> third point below
    assert _vmic_window(1.1, 0.5, 3) == [0.5, 1.0, 1.5]


def test_vmic_window_clamped_at_zero_extends_up():
    assert _vmic_window(0.2, 0.5, 4) == [0.0, 0.5, 1.0, 1.5]


@pytest.fixture
def atmos_grid():
    """A 2x2x1 (teff, logg) grid at feh=0, one atmosphere per node."""
    return {
        (5000.0, 4.0, 0.0): "t5000g40.mod",
        (6000.0, 4.0, 0.0): "t6000g40.mod",
        (5000.0, 4.5, 0.0): "t5000g45.mod",
        (6000.0, 4.5, 0.0): "t6000g45.mod",
    }


def test_select_nodes_representative_atmosphere_by_vmic():
    """4-tuple keys: per spatial node, the atmosphere nearest the query vmic wins."""
    grid = {
        (5000.0, 4.0, 0.0, 1.0): "v1.mod",
        (5000.0, 4.0, 0.0, 2.0): "v2.mod",
    }
    nodes, _ = select_nodes((5000.0, 4.0, 0.0, 1.9), grid)
    assert all(n[6] == "v2.mod" for n in nodes)
    nodes, _ = select_nodes((5000.0, 4.0, 0.0, 1.2), grid)
    assert all(n[6] == "v1.mod" for n in nodes)


def test_select_nodes_vmic_points_pchip_mode(atmos_grid):
    """vmic_points>2: n lattice levels per spatial node, node weight = spatial only."""
    nodes, vmic = select_nodes((5400.0, 4.0, 0.0, 1.3), atmos_grid, vmic_points=4)
    assert vmic == 1.3
    assert len(nodes) == 8  # 2 teff nodes x 4 vmic levels
    per_spatial = {}
    for jkey, teff, logg, feh, v, weight, path, absmet in nodes:
        per_spatial.setdefault(jkey[0], []).append((v, weight))
    for spatial, levels in per_spatial.items():
        vs = sorted(v for v, _ in levels)
        assert vs == [0.5, 1.0, 1.5, 2.0]
        ws = {w for _, w in levels}
        assert len(ws) == 1  # same (spatial) weight on every level
    # the spatial weights themselves still sum to 1 across nodes
    spatial_ws = [levels[0][1] for levels in per_spatial.values()]
    np.testing.assert_allclose(sum(spatial_ws), 1.0, rtol=1e-12)


def test_select_nodes_absmet_lookup(atmos_grid):
    nodes, _ = select_nodes((5000.0, 4.0, 0.0, 1.0), atmos_grid,
                            absmet_files={0.0: "absmet_m00"})
    assert nodes[0][7] == "absmet_m00"
    nodes, _ = select_nodes((5000.0, 4.0, 0.0, 1.0), atmos_grid,
                            absmet_files={-1.0: "absmet_m10"})
    assert nodes[0][7] == ""
