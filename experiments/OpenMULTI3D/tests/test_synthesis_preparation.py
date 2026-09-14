"""Common preparation and producing-run identity across synthesis entry points."""
from pathlib import Path
from types import SimpleNamespace

import pytest

from tumult.synthesis import SynthRequest, workflows as synth


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


def test_same_settings_for_single_and_batched_variants(calls):
    prefs = {'R': 100000., 'max_threads': 1}
    rv = {'mode': 'NLTE', 'abundances': 'Fe=0.2000', 'relative_abundances': True,
          'atoms': ['Fe.atom']}
    target, grid = (5750.,4.5,0.,1.), {'folder': 'grid', 'format': 'Marcs'}
    synth.compute(target,prefs,rv,5000.,5001.,grid=grid)
    synth.compute_batch(target,prefs,[rv],5000.,5001.,grid=grid)
    synth.compute_window_batch(target,prefs,[rv],[(5000.,5001.)],grid=grid)
    assert calls[0]['settings'] == calls[1]['settings_list'][0] == calls[2]['settings_list'][0]
    for key in ('launch_kwargs','star','atmos_files','nlte','vmic_step'):
        assert calls[0][key] == calls[1][key] == calls[2][key]


def test_batch_rejects_mixed_modes_and_negative_resolution(calls):
    args = ((5750.,4.5,0.,1.), {'max_threads': 1})
    with pytest.raises(ValueError, match='mode'):
        synth.compute_batch(*args,[{'mode':'LTE'},{'mode':'NLTE'}],grid={'folder':'grid'})
    with pytest.raises(ValueError, match='positive'):
        synth.compute_window_batch(args[0],{'R':-1.},[{}],[(5000.,5001.)],grid={'folder':'grid'})
    assert not calls


def test_request_owns_nested_inputs(calls):
    prefs, rv = {'max_threads':1}, {'atoms':['Fe.atom'],'teff':1.}
    request = SynthRequest(prefs,{'folder':'grid'},rv,(5750.,4.5,0.,1.))
    rv['atoms'].append('Ca.atom')
    prefs['R'] = -1.
    assert request.run_vals['atoms'] == ['Fe.atom']
    assert request.run_vals['teff'] == 5750.
    request.with_masks([{'fit': (5000., 5001.), 'synth': (5000., 5001.)}]).execute()
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


def test_uncropped_runs_keep_their_cache_identity(tmp_path, monkeypatch):
    """A per-node line-list crop must not move the folders of runs that have none.

    The digest only takes the crop into account when there is one, so what this pins is the
    condition for that: every way of asking for no crop returns None, and once there is one,
    its identity separates it from any other crop. Without the first, every cached run in
    every existing output directory becomes unreachable; without the second, two runs on one
    node with the same namelist and different line lists would share a folder.
    """
    from tumult.synthesis import api

    base = {'linelist_params': {'linelist_file': ['master.txt']},
            'spectrum_params': {'R': 300000.}}
    windows = [(5000., 5002.)]
    crop = {'skirt': 5.0, 'abundances': 'Fe=+0.40', 'relative_abundances': True}
    nodes = [(('n', 0), 5750., 4.5, 0.0, 1.0, 1.0, tmp_path / 'a.mod', None)]

    assert api._node_linelists('top', base, nodes, windows, None, {}) is None
    assert api._node_linelists('top', base, nodes, windows, {'skirt': 0.0}, {}) is None
    # a line-list folder cannot be cropped into one file
    folder = {**base, 'linelist_params': {'linelist_folder': 'lists/'}}
    assert api._node_linelists('top', folder, nodes, windows, crop, {}) is None

    written = []
    monkeypatch.setattr('tumult.prune.prune_batch',
                        lambda jobs, **kw: written.extend(Path(p) for p, _n in jobs) or written)
    mapping = api._node_linelists(tmp_path, base, nodes, windows, crop, {})
    assert list(mapping) == [('n', 0)]
    assert written == [Path(mapping[('n', 0)])]
    assert Path(mapping[('n', 0)]).name == 't5750g45m+00v10.txt'
    # the skirt, the windows and the rung the crop was made at all key the crop
    for changed in ({**crop, 'skirt': 2.0}, {**crop, 'abundances': 'Fe=+0.80'}):
        assert api._crop_digest(base, windows, changed) != api._crop_digest(base, windows, crop)
    assert api._crop_digest(base, [(5000., 5003.)], crop) != api._crop_digest(base, windows, crop)


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
