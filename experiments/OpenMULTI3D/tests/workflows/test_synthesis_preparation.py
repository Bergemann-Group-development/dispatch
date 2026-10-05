"""Common preparation and producing-run identity across synthesis entry points."""
from pathlib import Path
from types import SimpleNamespace

import pytest

from tumult.synthesis import AtmosphereGrid, SynthConfig, SynthRequest, workflows as synth


@pytest.fixture
def calls(monkeypatch, tmp_path):
    recorded = []
    monkeypatch.setattr(synth, 'discover_atmos_grid', lambda *a: {(5750,4.5,0.): tmp_path/'a.mod'})
    monkeypatch.setattr(synth, '_absmet_files', lambda *a: None)
    def record(**kwargs):
        recorded.append(kwargs)
        # execute_request unpacks (spectra per window, runs) from a window batch
        return ([[kwargs] * len(kwargs['windows'])], ()) if kwargs.get('return_runs') else kwargs
    for name in ('compute_spectrum', 'compute_spectra', '_compute_window_batch'):
        monkeypatch.setattr(synth, name, record)
    return recorded


def _request(tmp_path, **config):
    return SynthRequest(config=SynthConfig(root=tmp_path, **config),
                        atmosphere=AtmosphereGrid('grid'), star=(5750., 4.5, 0., 1.))


def test_same_settings_for_single_and_batched_variants(calls, tmp_path):
    request = _request(tmp_path, R_sampling=100000., threads=1)
    rv = {'mode': 'NLTE', 'abundances': 'Fe=0.2000', 'relative_abundances': True,
          'atoms': ['Fe.atom']}
    synth.compute(request, rv, 5000., 5001.)
    synth.compute_batch(request, [rv], 5000., 5001.)
    synth.compute_window_batch(request, [rv], [(5000., 5001.)])
    assert calls[0]['settings'] == calls[1]['settings_list'][0] == calls[2]['settings_list'][0]
    for key in ('launch_kwargs','star','atmos_files','nlte','vmic_step'):
        assert calls[0][key] == calls[1][key] == calls[2][key]


def test_batch_rejects_mixed_modes_and_negative_resolution(calls, tmp_path):
    request = _request(tmp_path, threads=1)
    with pytest.raises(ValueError, match='mode'):
        synth.compute_batch(request, [{'mode': 'LTE'}, {'mode': 'NLTE', 'atoms': ['a']}])
    with pytest.raises(ValueError, match='positive'):
        _request(tmp_path, R_sampling=-1.)
    assert not calls


def test_request_owns_nested_inputs(calls, tmp_path):
    preferences, rv = {'nx': 1}, {'atoms': ['Fe.atom'], 'mode': 'NLTE', 'teff': 1.}
    request = _request(tmp_path, threads=1, preferences=preferences)
    preferences['nx'] = 99
    assert request.config.preferences == {'nx': 1}
    resolved = synth.resolve_request(request.update(windows=[(5000., 5001.)]), rv)
    rv['atoms'].append('Ca.atom')
    assert resolved.values['atoms'] == ['Fe.atom']
    assert resolved.values['teff'] == 5750.
    synth.execute_request(resolved)
    assert len(calls) == 1 and calls[0]['windows'] == ((5000., 5001.),)


def test_cache_roots_follow_execution_directory(tmp_path, monkeypatch):
    execution, python_cwd = tmp_path/'dispatch', tmp_path/'client'
    execution.mkdir(); python_cwd.mkdir()
    monkeypatch.chdir(python_cwd)
    prefs = {'dispatch_cwd':str(execution),'synth_output_dir':'cache'}
    run = execution/'data/cache/bucket/t5750g45m+00v10'
    run.mkdir(parents=True)
    (run/'params.nml').write_text('')
    nml = execution/'cache/bucket/t5750g45m+00v10.nml'
    nml.parent.mkdir(parents=True); nml.write_text('')
    assert synth.list_cache(prefs)[0]['path'] == str(run)
    assert synth.delete_cache_entry(prefs,run)
    assert not run.exists() and not nml.exists()


def test_reference_linelist_uses_producing_run(tmp_path, monkeypatch):
    from tumult.linelist import LineList
    import pandas as pd
    own, other = tmp_path/'own', tmp_path/'newer'
    own.mkdir(); other.mkdir()
    (own/'linelist_out.txt').write_text('own')
    (other/'linelist_out.txt').write_text('unrelated newer dump')
    seen=[]
    def read(path):
        seen.append(Path(path))
        return pd.DataFrame({'ion':['Fe I'],'lam0':[5000.],'evi':[2.],'loggf':[-1.]})
    monkeypatch.setattr(LineList,'from_tumult_format',read)
    frame=synth.read_reference_linelist((own,))
    assert seen == [own/'linelist_out.txt']
    assert frame.lam0.iloc[0] == 5000.


def test_spectrum_table_preserves_first_intensity_and_units():
    import numpy as np
    import pandas as pd
    from astropy import units as u
    from tumult.korg import Korg
    from tumult.turbospectrum import Turbospectrum
    for cls in (Korg,Turbospectrum):
        result=cls.__new__(cls)
        result.df=pd.DataFrame([[5000.,.9,9.,.8,8.],[5001.,1.,10.,1.,10.]])
        assert result.get_spectrum(imu=0).values[1][0] == .8
        assert result.get_spectrum().values[1][0] == .9
        assert result.get_spectrum(norm=False).yunit == u.erg/(u.s*u.cm**2*u.AA)


def test_unpruned_runs_keep_their_cache_identity(tmp_path, monkeypatch):
    """Per-node line-list pruning must not move the folders of runs that have none.

    The digest takes the node lists into account only when they are a *selection*, a prune
    or a prepared list. Without that, either every cached run in every existing output
    directory becomes unreachable, or two runs on one node with the same namelist and
    different selections share a folder.
    """
    from tumult.synthesis import api

    base = {'linelist_params': {'linelist_file': ['master.txt']},
            'spectrum_params': {'R_sampling': 300000.}}
    windows = [(5000., 5002.)]
    prune = {'skirt': 5.0, 'abundances': 'Fe=+0.40', 'relative_abundances': True}
    nodes = [(('n', 0), 5750., 4.5, 0.0, 1.0, 1.0, tmp_path / 'a.mod', None)]

    # nothing named, nothing on disk: the readers stay on what the settings say
    assert api._node_linelists('top', base, nodes, windows, None, {}) is None
    folder = {**base, 'linelist_params': {'linelist_folder': 'lists/'}}
    assert api._node_linelists('top', folder, nodes, windows, prune, {}) is None

    written = []
    monkeypatch.setattr('tumult.prune.prune_batch',
                        lambda jobs, **kw: written.extend(Path(p) for p, _n in jobs) or written)
    mapping = api._node_linelists(tmp_path, base, nodes, windows, prune, {})
    assert list(mapping) == [('n', 0)]
    assert written == [Path(f) for f in mapping[('n', 0)]]
    assert [Path(f).name for f in mapping[('n', 0)]] == ['t5750g45m+00v10.txt']
    # one selection per node, at the composition the caller names; the prune reads the
    # master lists themselves, and hydrogen goes beside the selection whole
    written.clear()
    named = []
    monkeypatch.setattr('tumult.prune.prune_namelist',
                        lambda path, *a, **kw: named.append((Path(path).name, a, kw)) or {})
    master = tmp_path / 'master.txt'
    master.write_text("'Fe I'\n'0' 1\n 5000.0 1.0 -1.0\n'H I'\n'0' 1\n 6562.8 10.2 0.7\n"
                      "'Mg I'\n'0' 1\n 5172.7 2.7 -0.4\n")
    withh = {**base, 'linelist_params': {'linelist_file': [str(master)]}}
    mapping = api._node_linelists(tmp_path / 'h', withh, nodes, windows, prune, {})
    files = [Path(f) for f in mapping[('n', 0)]]
    assert [f.name for f in files] == ['hydrogen-0.txt', 't5750g45m+00v10.txt']
    assert files[0].read_text() == "'H I'\n'0' 1\n 6562.8 10.2 0.7\n"
    assert [(n, kw['abundances'], kw['exclude_elements']) for n, _a, kw in named] == [
        ('t5750g45m+00v10.txt', 'Fe=+0.40', 'H')]
    # the prune reads the master itself: HOP bisects it for the windows it is given
    sources, = [a[1] for _n, a, _kw in named]
    assert [Path(f) for f in sources] == [master]
    # the skirt, the windows and the rung the pruning was made at all key the pruning
    for changed in ({**prune, 'skirt': 2.0}, {**prune, 'abundances': 'Fe=+0.80'}):
        assert api._prune_digest(base, windows, changed) != api._prune_digest(base, windows, prune)
    assert api._prune_digest(base, [(5000., 5003.)], prune) != api._prune_digest(base, windows, prune)


def test_unfolded_ladder_is_one_batch_of_per_variant_runs(tmp_path, monkeypatch):
    """``fold=False`` must keep each variant's run identity and still schedule once.

    The point of not folding is that a rung's cache key is one abundance rather than a whole
    ladder, so a later solver iteration can reuse it; the point of batching is that the pool
    sees every (rung × node) run at once instead of one rung's nodes at a time. Both have to
    hold together, which is what this pins: distinct output folders per variant, one
    ``_run_batch`` call, and the folded call still issuing a single run per node.
    """
    from tumult.synthesis import api

    node = (('n', 0), 5750., 4.5, 0.0, 1.0, 1.0, str(tmp_path / 'a.mod'), None)
    monkeypatch.setattr(api, 'select_nodes', lambda *a, **kw: ([node], 1.0))
    monkeypatch.setattr(api, '_assemble_grid', lambda *a, **kw: object())
    monkeypatch.setattr(api, '_interp_crop', lambda grid, lo, hi, method: (lo, hi))
    batches = []

    def fake_run_batch(jobs, *a, **kw):
        batches.append([(key, str(path)) for key, path, _nml in jobs])
        return {key: SimpleNamespace(has_spectrum=True, mode='LTE', folder=path)
                for key, path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', fake_run_batch)
    variants = [{'composition_params': {'abundances': f'Fe={off:+.2f}'},
                 'linelist_params': {'linelist_file': ['master.txt']}}
                for off in (-0.25, 0.0, 0.25)]
    kwargs = dict(settings_list=variants, windows=[(5000., 5002.), (5100., 5102.)],
                  R=300000., launch_kwargs={'cwd': str(tmp_path)}, max_threads=1)

    api.compute_window_batch(tmp_path / 'synth', (5750., 4.5, 0., 1.), {}, fold=False,
                             **kwargs)
    assert len(batches) == 1, "the rungs must go out together, not one batch each"
    folders = [path for _key, path in batches[0]]
    assert len(folders) == len(variants) == len(set(folders))

    batches.clear()
    api.compute_window_batch(tmp_path / 'synth', (5750., 4.5, 0., 1.), {}, fold=True,
                             **kwargs)
    assert len(batches) == 1 and len(batches[0]) == 1, "a folded ladder is one run per node"


def test_node_store_launches_only_new_rungs_and_rejects_a_second_run(tmp_path, monkeypatch):
    """A workflow's node store must serve a rung it holds and refuse a conflicting run.

    A solver asks for overlapping ladders at every trial star; the rungs it already has on a
    node must come from memory, not from another launch or disk read. One grid point has one
    run in a workflow, so a request that would give the same point another run folder (here,
    other windows) is an error rather than a silent second entry.
    """
    from tumult.synthesis import api

    node = (('n', 0), 5750., 4.5, 0.0, 1.0, 1.0, str(tmp_path / 'a.mod'), None)
    monkeypatch.setattr(api, 'select_nodes', lambda *a, **kw: ([node], 1.0))
    monkeypatch.setattr(api, '_assemble_grid', lambda *a, **kw: object())
    monkeypatch.setattr(api, '_interp_crop', lambda grid, lo, hi, method: (lo, hi))
    launched = []

    def fake_run_batch(jobs, *a, **kw):
        launched.extend(nml['composition_params']['abundances'] for _key, _path, nml in jobs)
        return {key: SimpleNamespace(has_spectrum=True, mode='LTE', folder=path)
                for key, path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', fake_run_batch)
    store = api.NodeStore()

    def ladder(offsets, windows):
        variants = [{'composition_params': {'abundances': f'Fe={off:.4f}'},
                     'linelist_params': {'linelist_file': ['master.txt']}}
                    for off in offsets]
        return api.compute_window_batch(
            tmp_path / 'synth', (5750., 4.5, 0., 1.), {}, settings_list=variants,
            windows=windows, R=300000., launch_kwargs={'cwd': str(tmp_path)},
            max_threads=1, store=store, abundances=list(offsets))

    windows = [(5000., 5002.), (5100., 5102.)]
    ladder((-0.25, 0.0, 0.25), windows)
    assert len(launched) == 3
    launched.clear()
    results = ladder((0.0, 0.25, 0.5), windows)
    assert launched == ['Fe=0.5'], "only the new rung is launched"
    assert len(results) == 3 and all(len(r) == len(windows) for r in results)
    with pytest.raises(RuntimeError, match="Two runs for one grid point"):
        ladder((0.0,), [(5000., 5002.)])


def test_probe_stars_share_the_target_runs_and_add_only_their_own_nodes(tmp_path, monkeypatch):
    """A probe star is interpolated from the same variants' runs; a node only it needs is
    launched with the target's, once, and every star gets its own interpolation."""
    from tumult.synthesis import api

    shared = (('a', 1.0), 5750., 4.5, 0.0, 1.0, 1.0, str(tmp_path / 'a.mod'), None)
    extra = (('b', 1.0), 6000., 4.5, 0.0, 1.0, 1.0, str(tmp_path / 'b.mod'), None)

    def select_nodes(star, *a, **kw):
        return ([shared] if star[0] < 5800. else [shared, extra]), star[3]

    monkeypatch.setattr(api, 'select_nodes', select_nodes)
    monkeypatch.setattr(api, '_assemble_grid',
                        lambda views, star, *a, **kw: (star[0], sorted(views)))
    monkeypatch.setattr(api, '_interp_crop', lambda grid, lo, hi, method: grid)
    launched = []

    def fake_run_batch(jobs, *a, **kw):
        launched.extend(key for key, _path, _nml in jobs)
        return {key: SimpleNamespace(has_spectrum=True, mode='LTE', folder=path)
                for key, path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', fake_run_batch)
    variants = [{'composition_params': {'abundances': f'Fe={off:.4f}'},
                 'linelist_params': {'linelist_file': ['master.txt']}} for off in (0.0, 0.5)]
    results = api.compute_window_batch(
        tmp_path / 'synth', (5750., 4.5, 0., 1.), {}, settings_list=variants,
        windows=[(5000., 5002.)], R=300000., launch_kwargs={'cwd': str(tmp_path)},
        max_threads=1, store=api.NodeStore(), abundances=[0.0, 0.5],
        probes=[(5850., 4.5, 0., 1.)])
    assert len(launched) == 4, "two variants on the union of two nodes, each once"
    target, probe = results
    assert target[0][0] == (5750., [('a', 1.0)])
    assert probe[1][0] == (5850., [('a', 1.0), ('b', 1.0)])
