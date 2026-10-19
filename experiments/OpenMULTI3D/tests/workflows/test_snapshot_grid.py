"""Registering and interpolating a 3D snapshot grid through the generic path pattern."""
from pathlib import Path
import re

import pytest

from tumult.synthesis import snapshot_grid, workflows
from tumult.synthesis.atmos_select import select_nodes


def _node(folder, name, snaps=(8, 31, 32), model=None):
    model = model or name.replace('_new', '')
    node = folder / name
    node.mkdir(parents=True)
    for snap in snaps:
        (node / f'atm3d.{model}05_{snap:05d}_80x80x240').write_text('')
        (node / f'mesh.{model}05_{snap:05d}_80x80x240').write_text('')
    (node / 'EOSrhoe.tab').write_text('')
    return node


@pytest.fixture
def grid(tmp_path):
    for name in ('t50g45m00', 't55g45m00', 't50g45m10', 't55g45m10_new', 't5777g44m30'):
        _node(tmp_path, name)
    (tmp_path / 'README').write_text('')
    return tmp_path


def _snap(name):
    return re.search(r'_(\d{5})_', name).group(1)


def _discover(folder, snapshots='1'):
    return workflows.discover_atmos_grid(folder, workflows.STAGGER_PATTERN, snapshots=snapshots)


def test_pattern_grammar_covers_the_stagger_conventions():
    p = workflows.compile_pattern(workflows.STAGGER_PATTERN)
    assert p.search('t50g20m00/atm3d.t50g20m0005_00008_80x80x240') == {
        'Teff': 5000., 'logg': 2.0, 'FeH': 0.0, 'snap': 8.}
    assert p.search('t50g47p05_new/atm3d.t50g47p0505_00031_80x80x240') == {
        'Teff': 5000., 'logg': 4.7, 'FeH': 0.5, 'snap': 31.}
    assert p.search('t50g45m00/mesh.t50g45m0005_00008_80x80x240') is None
    assert p.search('t50g45m00/EOSrhoe.tab') is None
    # the MARCS default is unchanged by the extended grammar
    assert workflows.recover_params(
        workflows.DEFAULT_PATTERN, 'p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00.mod') == (
        5777., 4.4, 0.0, 1.0)


def test_pattern_format_items_constrain_the_match():
    # the format that wrote the name decides what a token may consume
    two = workflows.compile_pattern('t{teff/100:d}g{logg*10:.0f}{feh*10:+03.0f}')
    assert two.search('atm3d.t50g45m0005_00008') == {'Teff': 5000., 'logg': 4.5, 'FeH': 0.0}
    loose = workflows.compile_pattern('t{teff/100:d}g{logg*10:.0f}{feh*10:+.0f}')
    assert loose.search('atm3d.t50g45m0005_00008')['FeH'] == -0.5   # the version digits leak in
    signed = workflows.compile_pattern('z{FeH:+.2f}_g{logg}_t{Teff}')
    assert signed.search('z-1.50_g4.5_t5000')['FeH'] == -1.5
    assert signed.search('zp1.50_g4.5_t5000')['FeH'] == 1.5
    assert signed.search('z1.50_g4.5_t5000') is None                 # sign is mandatory
    with pytest.raises(ValueError, match='different fields'):
        workflows.compile_pattern('t{Teff}g{logg}m{FeH}v{vmic}|t{Teff}g{logg}m{FeH}')
    with pytest.raises(ValueError, match='Unrecognised format item'):
        workflows.compile_pattern('t{Teff:digits=2}g{logg}m{FeH}')


def test_recursive_discovery_selects_snapshots_per_node_and_prefers_new_models(grid):
    _node(grid, 't55g45m10', snaps=(1,))            # superseded by t55g45m10_new
    nodes = _discover(grid, '2')
    assert set(nodes) == {(5000., 4.5, 0.), (5500., 4.5, 0.), (5000., 4.5, -1.),
                          (5500., 4.5, -1.), (5777., 4.4, -3.)}
    # a count picks the leading entries of snapshot_grid.halton_order, in snapshot order
    assert [p.name for p in nodes[(5500., 4.5, -1.)]] == [
        'atm3d.t55g45m1005_00008_80x80x240', 'atm3d.t55g45m1005_00032_80x80x240']
    assert len(_discover(grid, 'all')[(5000., 4.5, 0.)]) == 3
    picked = _discover(grid, '8, 32')[(5000., 4.5, 0.)]
    assert [_snap(p.name) for p in picked] == ['00008', '00032']
    with pytest.raises(ValueError, match='no snapshots'):
        _discover(grid, '8,99')


def test_flat_grid_discovery_keeps_one_file_per_node(tmp_path):
    for name in ('p5000_g+4.5_m0.0_t01_st_z+0.00_a+0.00.mod',
                 'p5000_g+4.5_m0.0_t01_st_z+0.00_a+0.40.mod',
                 'p5500_g+4.5_m0.0_t01_st_z+0.00_a+0.00.mod'):
        (tmp_path / name).write_text('')
    files = workflows.discover_atmos_grid(tmp_path)
    assert not snapshot_grid.is_snapshot_node(files)
    assert files[(5000., 4.5, 0.0, 1.0)].name.endswith('a+0.40.mod')   # last wins


def test_snapshot_nodes_have_no_vmic_lattice_and_share_the_spatial_weight(grid):
    files = _discover(grid, '2')
    nodes, vmic = select_nodes((5250., 4.5, 0.0, 1.3), files, method='Delaunay',
                               vmic_step=0.5)
    assert vmic == 1.3
    assert all(node[4] is None for node in nodes)
    assert abs(sum(node[5] for node in nodes) - 1.0) < 1e-12
    spatial = {node[0][0] for node in nodes}
    assert spatial == {(5000., 4.5, 0.), (5500., 4.5, 0.)}
    assert {_snap(node[0][1]) for node in nodes} == {'00008', '00032'}
    assert all(Path(node[6]).name.startswith('atm3d.') for node in nodes)


def test_snapshot_node_jobs_name_each_snapshot_and_write_no_vmic(grid):
    from tumult.synthesis import api
    files = _discover(grid, '2')
    nodes, _ = select_nodes((5000., 4.5, 0.0, 1.0), files)
    jobs, keys, weights, params = api._spec_jobs({'atmos_params': {'atmos_format': 'Stagger2'}},
                                                 nodes, Path('synth'), False, {})
    names = [path.name for _key, path, _nml in jobs]
    assert names == ['t5000g450m+00_atm3d_t50g45m0005_00008_80x80x240',
                     't5000g450m+00_atm3d_t50g45m0005_00032_80x80x240']
    for _key, _path, nml in jobs:
        assert 'vmic' not in nml['atmos_params']
        assert nml['atmos_params']['atmos_format'] == 'Stagger2'
    assert weights == [0.5, 0.5] and all(p[3] == 0.0 for p in params)
    dep = api._dep_jobs({}, nodes, Path('dep'))
    assert [path.name for _k, path, _n in dep] == names
    assert all('vmic' not in nml['atmos_params'] for _k, _p, nml in dep)


def test_grid_record_discovery_and_mesh_route_by_format(grid, tmp_path):
    from tumult.synthesis import workflows
    record = {'folder': str(grid), 'format': 'Stagger2', 'snapshots': '1',
              'pattern': workflows.STAGGER_PATTERN, 'resolution': [80, 80, 240]}
    files = snapshot_grid.discover_grid(record)
    assert snapshot_grid.is_snapshot_node(files)
    assert workflows._grid_mesh(record) == {'nx': 80, 'ny': 80, 'nz': 240}
    assert workflows._grid_mesh({'folder': str(grid), 'format': 'Marcs'}) == {}
    with pytest.raises(ValueError):
        workflows._grid_mesh({'resolution': [80, 80]})


def test_prepare_builds_snapshot_grid_settings(grid, tmp_path):
    from tumult.synthesis import workflows
    from tumult.synthesis import AtmosphereGrid, SynthConfig, SynthRequest
    atmosphere = AtmosphereGrid(grid, pattern=workflows.STAGGER_PATTERN, format='Stagger2',
                                snapshots='1', resolution=(40, 40, 120))
    request = SynthRequest(config=SynthConfig(root=tmp_path, R_sampling=50000),
                           atmosphere=atmosphere, star=(5250., 4.5, 0.0, 1.0))
    common, settings, _ = workflows._prepare(request, [{}])
    assert snapshot_grid.is_snapshot_node(common['atmos_files'])
    assert settings[0]['atmos_params'] == {'atmos_format': 'Stagger2',
                                           'nx': 40, 'ny': 40, 'nz': 120}


def test_grid_assembly_averages_a_node_over_its_snapshots(grid):
    """The interpolation grid sees one spectrum per stellar point, the snapshot mean."""
    import numpy as np
    from tumult import Spectrum
    from tumult.synthesis import api
    files = _discover(grid, '2')
    nodes, vmic = select_nodes((5250., 4.5, 0.0, 1.0), files)
    keys = [node[0] for node in nodes]

    class Run:
        has_spectrum, mode, folder = True, 'LTE', grid

        def __init__(self, level):
            self.level = level

        def crop(self, **_kwargs):
            return Spectrum(np.array([5250., 5251.]), np.full(2, self.level))

    # node 5000 K: snapshots at 0.2 and 0.4 (mean 0.3); node 5500 K: 0.8 and 1.0 (mean 0.9)
    runs = {key: Run(0.2 + 0.2 * i if key[0][0] == 5000. else 0.8 + 0.2 * i)
            for i, key in enumerate(k for k in keys if k[0][0] == 5000.)}
    runs.update({key: Run(0.8 + 0.2 * i)
                 for i, key in enumerate(k for k in keys if k[0][0] == 5500.)})
    weights = [node[5] for node in nodes]
    params = [(n[1], n[2], n[3], 0.0) for n in nodes]
    g = api._assemble_grid(runs, (5250., 4.5, 0.0), vmic, keys, weights, params, 2)
    assert len(g.spectra) == 2 and np.allclose(g.interp_weights, [0.5, 0.5])
    spec = api._interp_crop(g, 5250., 5251., 'Delaunay')
    assert np.allclose(spec.values[1], 0.6)
