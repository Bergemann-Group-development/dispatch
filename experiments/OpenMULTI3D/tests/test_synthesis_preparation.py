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
        return kwargs
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
    request.compute(5000.,5001.)
    assert len(calls) == 1


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
