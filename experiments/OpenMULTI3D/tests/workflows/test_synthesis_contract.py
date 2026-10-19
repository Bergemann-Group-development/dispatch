"""Public request, preparation and departure contracts.

The durable-submission contracts, which need the workflow daemon, are in
``tests/execution/test_synthesis_submission.py``.
"""
from pathlib import Path

import numpy as np
import pytest

from tumult import Spectrum
from tumult.synthesis import (Atmosphere, AtmosphereGrid, ResolvedSynthRequest, Snapshot,
                              SynthConfig, SynthesisResult, SynthRequest, workflows)


def _request(atmosphere, tmp_path, windows=(), star=(5750, 4.5, 0, 1), **config):
    config = SynthConfig(**{'root': tmp_path, 'cache': 'out', 'R_sampling': 50000,
                            **config})
    return SynthRequest(config=config, atmosphere=atmosphere, star=star, windows=windows)


def test_all_atmosphere_kinds_roundtrip(tmp_path):
    requests = [_request(AtmosphereGrid('grid'), tmp_path, [(4999.5, 5001.5)]),
                _request(Atmosphere('one.mod'), tmp_path, [(4999.5, 5001.5)]),
                _request(Snapshot('cube.dat', 'Stagger', (8, 8, 16)), tmp_path,
                         [(4999.5, 5001.5)], star=(5750, 4.5, 0))]
    for request in requests:
        reopened = SynthRequest.from_json(request.to_json())
        assert type(reopened.atmosphere) is type(request.atmosphere)
        assert reopened.windows == ((4999.5, 5001.5),)
        assert reopened.to_dict() == request.to_dict()


def test_resolve_uses_request_windows_and_maintained_dispatch_padding(tmp_path):
    from tumult.abundance import SYNTH_PAD_AA
    request = _request(Atmosphere('one.mod'), tmp_path, [(4999.5, 5001.5)])
    resolved = workflows.resolve_request(request, {'only_elements': 'Mg'})
    assert resolved.windows == ((4999.5, 5001.5),)
    assert resolved.padded_windows == ((4999.5 - SYNTH_PAD_AA, 5001.5 + SYNTH_PAD_AA),)
    assert Path(resolved.execution_root) == tmp_path.resolve()
    assert resolved.values['only_elements'] == 'Mg' and resolved.values['teff'] == 5750
    assert ResolvedSynthRequest.from_json(resolved.to_json()).to_dict() == resolved.to_dict()


def test_preparation_cache_tracks_changed_inputs(tmp_path, monkeypatch):
    from tumult import linestrength
    master, atmos = tmp_path / 'master.txt', tmp_path / 'a.mod'
    master.write_text('one')
    atmos.write_text('atmos')
    calls = []

    def prune_batch(jobs, **_kwargs):
        for path, _nml in jobs:
            calls.append(Path(path))
            Path(path).write_text('prepared')
        return [path for path, _nml in jobs]

    monkeypatch.setattr('tumult.prune.prune_batch', prune_batch)
    a = linestrength.prepare_linelist(tmp_path / 'cache', [atmos], master,
                                      [(5000, 5002)], R=50000)
    b = linestrength.prepare_linelist(tmp_path / 'cache', [atmos], master,
                                      [(5000, 5002)], R=50000)
    assert a.folder == b.folder and len(calls) == 1
    master.write_text('changed master input')
    c = linestrength.prepare_linelist(tmp_path / 'cache', [atmos], master,
                                      [(5000, 5002)], R=50000)
    assert c.folder != a.folder and len(calls) == 2


def test_run_settings_name_a_linelist_folders_files_explicitly(tmp_path):
    """A line-list folder reaches dispatch as its sorted files, never as
    ``linelist_folder``: dispatch lists a folder through a shell temp file in the shared
    working directory, and the concurrent runs of one synthesis delete it from under each
    other (``Cannot open file '.temp_files000'``)."""
    from tumult.synthesis import workflows
    folder = tmp_path / 'lists'
    folder.mkdir()
    (folder / 'Mg.txt').write_text("'Mg I'\n")
    (folder / 'Fe.txt').write_text("'Fe I'\n")
    (folder / '.ignored').write_text('metadata')
    assert workflows._linelist_settings(folder) == {
        'linelist_file': [str(folder / 'Fe.txt'), str(folder / 'Mg.txt')]}
    assert workflows._linelist_settings(folder / 'Fe.txt') == {
        'linelist_file': [str(folder / 'Fe.txt')]}
    assert workflows._linelist_settings('') == {}


def test_preparation_expands_a_linelist_folder_for_the_standalone_pruner(tmp_path,
                                                                        monkeypatch):
    from tumult import linestrength
    master, atmos = tmp_path / 'master', tmp_path / 'a.mod'
    master.mkdir()
    (master / 'Fe.txt').write_text("'Fe I'\n")
    (master / 'Mg.txt').write_text("'Mg I'\n")
    (master / '.ignored').write_text('metadata')
    atmos.write_text('atmos')
    seen = []

    def prune_batch(jobs, **_kwargs):
        for path, nml in jobs:
            seen.append(nml['linelist_params']['linelist_file'])
            Path(path).write_text('prepared')

    monkeypatch.setattr('tumult.prune.prune_batch', prune_batch)
    linestrength.prepare_linelist(tmp_path / 'cache', [atmos], master,
                                  [(5000, 5002)], R=50000)
    assert seen == [[str(master / 'Fe.txt'), str(master / 'Mg.txt')]]


def test_individual_atmosphere_resolves_nearest_marcs_reference(tmp_path, monkeypatch):
    from tumult.synthesis import workflows
    marcs = tmp_path / 'marcs'
    marcs.mkdir()
    cool, solar = marcs / 'cool.mod', marcs / 'solar.mod'
    cool.write_text('cool')
    solar.write_text('solar')
    monkeypatch.setattr(workflows, 'discover_atmos_grid', lambda *args: {
        (4500., 2.0, -1.0): cool, (5750., 4.5, 0.0): solar})
    atmosphere = Atmosphere('science.mod', 'Text', preparation=AtmosphereGrid(marcs))
    resolved = workflows.resolve_request(
        _request(atmosphere, tmp_path, star=(5700, 4.4, 0.0, 1.0)))
    assert resolved.preparation_references == (str(solar.resolve()),)
    assert resolved.preparation_aliases[str((tmp_path / 'science.mod').resolve())] == \
        str(solar.resolve())


def test_snapshot_execution_routes_directly_without_grid_or_vmic(tmp_path, monkeypatch):
    """A snapshot is one run over every window, read through a window file."""
    from tumult.synthesis import snapshot, workflows
    run = type('Run', (), {'folder': tmp_path / 'run',
                           'crop': lambda self, **kwargs: ('spectrum', kwargs)})()
    seen = []

    def direct(top, path, **kwargs):
        seen.append((top, Path(path), kwargs))
        return run

    monkeypatch.setattr(snapshot, 'synthesise_snapshot', direct)
    windows = [(4999.5, 5001.5), (6100.0, 6102.0)]
    request = _request(Snapshot('cube.dat', 'Stagger', (4, 5, 6)), tmp_path,
                       windows, star=(5750, 4.5, 0))
    result = workflows.execute_request(request)
    assert result.geometry == '3d' and len(seen) == 1
    assert [kwargs['window'] for _name, kwargs in result.spectra] == windows
    assert seen[0][1] == (tmp_path / 'cube.dat').resolve()
    assert seen[0][2]['resolution'] == (4, 5, 6)
    assert 'vmic_step' not in seen[0][2]
    spectrum = seen[0][2]['settings']['spectrum_params']
    assert 'aa_blue' not in spectrum
    assert len(Path(spectrum['window_file']).read_text().splitlines()) == 3


def test_saved_departures_can_be_opened_from_native_artifacts(tmp_path):
    from tumult.synthesis import Departures
    folder = tmp_path / 'dep'
    (folder / 'save').mkdir(parents=True)
    (folder / 'save' / 'pop_Mg_meta.txt').write_text('')
    (folder / 'params.nml').write_text(
        "&atmos_params\n atmos_file='cube.dat' nx=4 ny=5 nz=6\n/\n"
        "&composition_params\n abundances='Mg=0.2'\n/\n")
    field = Departures.load(folder)
    assert field.saved and field.resolution == (4, 5, 6)
    assert field.species == ('Mg',) and field.composition['abundances'] == 'Mg=0.2'


def test_window_batch_reports_runs_without_a_ladder(tmp_path, monkeypatch):
    """A single variant without abundance overrides is not folded, yet still names its runs."""
    from tumult.synthesis import api

    class Run:
        def __init__(self, folder):
            self.folder = folder

    def fake_run_batch(jobs, *_args, **_kwargs):
        return {key: Run(tmp_path / str(path)) for key, path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', fake_run_batch)
    monkeypatch.setattr(api, '_assemble_grid', lambda *a, **k: 'grid')
    monkeypatch.setattr(api, '_interp_crop', lambda grid, lo, hi, method: (lo, hi))
    grid = tmp_path / 'grid'
    grid.mkdir()
    (grid / 'a.mod').write_text('')
    settings = [{'linelist_params': {'linelist_file': ['master.txt']}}]
    values, runs = api.compute_window_batch(
        tmp_path / 'synth', (5750., 4.5, 0., 1.), {(5750., 4.5, 0.): grid / 'a.mod'},
        settings_list=settings, windows=[(5000., 5002.), (5100., 5102.)], R=300000.,
        launch_kwargs={'cwd': str(tmp_path)}, max_threads=1, return_runs=True)
    assert values == [[(5000., 5002.), (5100., 5102.)]]
    assert len(runs) == 2 and len(set(runs)) == 2


# --- unsupported combinations are refused before submission ----------------------------
def _saved_departures(tmp_path):
    from tumult.synthesis import Departures
    folder = tmp_path / 'dep'
    (folder / 'save').mkdir(parents=True)
    (folder / 'save' / 'pop_Mg_meta.txt').write_text('')
    (folder / 'params.nml').write_text(
        "&atmos_params\n atmos_file='one.mod' nx=1 ny=1 nz=6\n/\n"
        "&atom_params\n atom_file='mg.atom'\n/\n"
        "&composition_params\n abundances='Mg=0.2'\n/\n")
    return Departures.load(folder)


def test_unsupported_combinations_are_refused_before_a_record_exists(tmp_path):
    from tumult.departures import DepartureRequest
    from tumult.synthesis import submission
    field = _saved_departures(tmp_path)
    one = _request(Atmosphere('one.mod'), tmp_path, [(5000, 5001)])
    supported = one.update(departures=field)
    assert workflows.resolve_request(supported).values['mode'] == 'NLTE'
    with pytest.raises(ValueError, match='DepartureGrid'):
        _request(AtmosphereGrid('grid'), tmp_path).update(departures=field)
    with pytest.raises(ValueError, match='model atoms'):
        workflows.resolve_request(one, {'mode': 'NLTE'})
    with pytest.raises(ValueError, match='prepared'):
        workflows.resolve_request(one.update(departures=DepartureRequest('mg.atom')))
    with pytest.raises(ValueError, match='unknown operation'):
        submission.submit_operation('synthesis.compute', one, root=tmp_path / 'jobs')
    assert not (tmp_path / 'jobs').exists()      # nothing was recorded for a refused request


def test_fit_refuses_what_fitting_does_not_support(tmp_path):
    from tumult.synthesis import check_fit
    windows = [(5000, 5001), (5002, 5003)]
    grid = _request(AtmosphereGrid('grid'), tmp_path, windows)
    check_fit(grid, fit_stellar=['teff'], abund_fits={'Mg': 0.1}, fit_post={'vshift': 0})

    one = _request(Atmosphere('one.mod'), tmp_path, windows)
    with pytest.raises(ValueError, match='1D atmospheres and grids'):
        check_fit(_request(Snapshot('cube.dat', 'Stagger', (4, 4, 8)), tmp_path, windows,
                           star=(5750, 4.5, 0)), fit_post={'vshift': 0})
    with pytest.raises(ValueError, match='Frozen departures'):
        check_fit(one.update(departures=_saved_departures(tmp_path)), abund_fits={'Mg': 0.0})
    with pytest.raises(ValueError, match='single atmosphere file'):
        check_fit(one, fit_stellar=['teff'])
    with pytest.raises(ValueError, match='at least one parameter'):
        check_fit(one)
    with pytest.raises(ValueError, match='window'):
        check_fit(one.update(windows=()), abund_fits={'Mg': 0.0})


def test_fit_synthesises_each_trial_point_as_the_request_over_the_span(tmp_path, monkeypatch):
    """Every trial point is the fit's own request at that point, over one span window, so
    the configuration, expert settings, sampling and padding reach dispatch unchanged; the
    observation is compared only inside the fit windows."""
    from tumult.synthesis import Abundances, api, operations, snapshot_grid
    request = _request(AtmosphereGrid('grid'), tmp_path, [(5000, 5001), (5002, 5003)],
                       abundance_step=0.5, preferences={'fit_max_rounds': 3},
                       expert={'spectrum_params': {'daa': 0.01}})
    request = request.update(abundances={'Ca': 0.1})
    executed = []

    def synthesize(trial, background_lines, cancel_token):
        executed.append(trial)
        xx = np.linspace(5000., 5003., 31)
        return SynthesisResult((Spectrum(xx, np.ones_like(xx)),), trial.windows, trial, {}, {})

    nodes = {(5500., 4.5, 0., 1.): 'a', (6000., 4.5, 0., 1.): 'b'}
    monkeypatch.setattr(operations, '_synthesize', synthesize)
    monkeypatch.setattr(operations, 'preflight', lambda request: None)
    monkeypatch.setattr(snapshot_grid, 'discover_grid', lambda grid, root=None: nodes)
    seen = {}

    def adaptive(target, synth_point_fn, axes, guesses, **kwargs):
        seen.update(target=target, axes=axes, guesses=guesses, **kwargs)
        for teff in (5500., 6000.):
            synth_point_fn({'teff': teff, 'Mg': 0.25})
        best = {'teff': 5600., 'Mg': 0.25, 'vshift': 1.5}
        return best, synth_point_fn(best), None

    monkeypatch.setattr(api, 'adaptive_fit', adaptive)
    xx = np.arange(4999., 5004., 0.1)
    observation = Spectrum(xx, np.ones_like(xx))
    result = operations.fit(request, target=observation,
                            fit_windows=[(5000.2, 5000.8), (5002.2, 5002.6)],
                            fit_stellar=['teff'], abund_fits={'Mg': 0.0},
                            fit_post={'vshift': 0.0})

    assert result.best['teff'] == 5600. and result.post == {'vshift': 1.5}
    assert seen['axes'] == ['teff', 'Mg'] and seen['post_fit'] == ['vshift']
    assert seen['axis_nodes'] == {'teff': [5500., 6000.]} and seen['spacings'] == {'Mg': 0.5}
    assert seen['max_rounds'] == 3
    inside = seen['target'].values[0]
    assert inside.min() >= 5000.2 and inside.max() <= 5002.6
    assert not ((inside > 5000.8) & (inside < 5002.2)).any()   # the gap is dropped
    # the trial points (and the raw recompute of the best one): the request at each point
    # over the one span window
    trial = executed[0]
    assert [v.star.teff for v in executed] == [5500., 6000., 5600., 5600.]
    assert trial.windows == ((5000., 5003.),)
    assert trial.abundances == Abundances({'Ca': 0.1, 'Mg': 0.25})
    assert operations.run_values(trial)['abundances'] == 'Ca=0.1000 Mg=0.2500'
    assert trial.config == request.config and trial.atmosphere == request.atmosphere
    spectrum = result.spectrum
    assert spectrum.values[0][0] == 5000. and spectrum.values[0][-1] == 5003.


@pytest.mark.parametrize('with_linelist', [False, True])
def test_snapshot_correction_uses_atom_grid_without_unbounded_R_sampling(
        tmp_path, with_linelist):
    from tumult.synthesis import snapshot

    atom = tmp_path / 'atom.ba'
    line_grid = tmp_path / 'line_grid.txt'
    request = _request(Snapshot('cube.atmos', 'M3DIS', (20, 20, 64)), tmp_path,
                       star=(4548, 3.964, 0), R_sampling=1e6, linelist_skirt=0)
    run_values = {'mode': 'NLTE', 'atoms': [str(atom)], 'atom_lines': True,
                  'only_atom': not with_linelist, 'line_grid': str(line_grid)}
    resolved = workflows.resolve_request(request, run_values)
    for section in ('rt_iteration_params', 'rt_output_params'):
        assert resolved.settings[section]['sample_back_opac'] is True
        assert resolved.settings[section]['sample_linelist'] is True
    field = snapshot.Departures(
        folder=tmp_path / 'departures', name=Path('synth/dep'),
        atmosphere='cube.atmos', resolution=(20, 20, 64),
        settings=snapshot.departure_settings(resolved.settings, 'cube.atmos'))
    formal = snapshot.spectrum_settings(field, resolved.settings)
    assert not formal.get('spectrum_params', {}).get('R_sampling')
    assert formal['atom_params']['line_grid_file'] == str(line_grid)
    assert formal['atmos_params']['read_departures']
    assert formal['rt_iteration_params']['maxiter'] == 0

    # A user-selected wavelength window still needs the requested sampling.
    windowed = workflows.resolve_request(request.update(windows=[(4552, 4556)]), run_values)
    assert windowed.settings['spectrum_params']['R_sampling'] == 1e6


def test_snapshot_sampled_opacity_defaults_allow_expert_override(tmp_path):
    expert = {
        'rt_iteration_params': {'sample_back_opac': False, 'sample_linelist': False},
        'rt_output_params': {'sample_back_opac': False, 'sample_linelist': False},
    }
    request = _request(Snapshot('cube.atmos', 'M3DIS', (20, 20, 64)), tmp_path,
                       star=(4548, 3.964, 0), expert=expert)
    settings = workflows.resolve_request(request, {}).settings
    for section in ('rt_iteration_params', 'rt_output_params'):
        assert settings[section]['sample_back_opac'] is False
        assert settings[section]['sample_linelist'] is False

    one_d = workflows.resolve_request(_request(Atmosphere('one.mod'), tmp_path), {}).settings
    assert 'rt_iteration_params' not in one_d
    assert 'rt_output_params' not in one_d
