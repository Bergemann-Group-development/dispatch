"""Public request, preparation, departure and durable-submission contracts."""
import time
from pathlib import Path

import numpy as np
import pytest

from tumult import Spectrum
from tumult.synthesis import (Atmosphere1D, GridTarget, MaskWindow, ResolvedSynthRequest, Snapshot3D,
                              SynthesisResult, NativeRequest)


def _request(selection, tmp_path, masks=()):
    return NativeRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out', 'R': 50000},
                        run_vals={'mode': 'LTE'}, selection=selection, masks=masks)


def test_all_atmosphere_selections_roundtrip_and_capture_masks(tmp_path):
    masks = [{'fit': [5000, 5001], 'synth': [4999.5, 5001.5], 'include': True}]
    selections = [GridTarget({'folder': 'grid'}, (5750, 4.5, 0, 1)),
                  Atmosphere1D('one.mod', (5750, 4.5, 0, 1)),
                  Snapshot3D('cube.dat', (5750, 4.5, 0), (8, 8, 16))]
    for selection in selections:
        request = _request(selection, tmp_path, masks)
        masks[0]['synth'][0] = 1
        reopened = NativeRequest.from_json(request.to_json())
        assert type(reopened.selection) is type(selection)
        assert reopened.masks[0].synth == (4999.5, 5001.5)
        masks[0]['synth'][0] = 4999.5


def test_resolve_uses_curated_windows_and_maintained_dispatch_padding(tmp_path):
    from tumult.abundance import SYNTH_PAD_AA
    request = _request(Atmosphere1D('one.mod', (5750, 4.5, 0, 1)), tmp_path,
                       [MaskWindow((5000, 5001), (4999.5, 5001.5))])
    resolved = request.resolve()
    assert resolved.windows == ((4999.5, 5001.5),)
    assert resolved.padded_windows == ((4999.5 - SYNTH_PAD_AA, 5001.5 + SYNTH_PAD_AA),)
    assert Path(resolved.execution_root) == tmp_path.resolve()
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


def test_individual_atmosphere_resolves_nearest_marcs_reference(tmp_path, monkeypatch):
    from tumult.synthesis import workflows
    marcs = tmp_path / 'marcs'
    marcs.mkdir()
    cool, solar = marcs / 'cool.mod', marcs / 'solar.mod'
    cool.write_text('cool')
    solar.write_text('solar')
    monkeypatch.setattr(workflows, 'discover_atmos_grid', lambda *args: {
        (4500., 2.0, -1.0): cool, (5750., 4.5, 0.0): solar})
    selection = Atmosphere1D('science.mod', (5700, 4.4, 0.0, 1.0), 'Text',
                             {'folder': str(marcs), 'format': 'Marcs'})
    resolved = _request(selection, tmp_path).resolve()
    assert resolved.preparation_references == (str(solar.resolve()),)
    assert resolved.preparation_aliases[str((tmp_path / 'science.mod').resolve())] == \
        str(solar.resolve())


def test_snapshot_execution_routes_directly_without_grid_or_vmic(tmp_path, monkeypatch):
    from tumult.synthesis import snapshot, workflows
    run = type('Run', (), {'folder': tmp_path / 'run',
                           'crop': lambda self, **kwargs: ('spectrum', kwargs)})()
    seen = []

    def direct(top, path, **kwargs):
        seen.append((top, Path(path), kwargs))
        return run

    monkeypatch.setattr(snapshot, 'synthesise_snapshot', direct)
    request = _request(Snapshot3D('cube.dat', (5750, 4.5, 0), (4, 5, 6)), tmp_path,
                       [MaskWindow((5000, 5001), (4999.5, 5001.5))])
    result = workflows.execute_request(request)
    assert result.geometry == '3d' and len(result.spectra) == 1
    assert seen[0][1] == (tmp_path / 'cube.dat').resolve()
    assert seen[0][2]['resolution'] == (4, 5, 6)
    assert 'vmic_step' not in seen[0][2]


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


def test_durable_submission_reopens_saved_result(tmp_path, monkeypatch, workflow_daemon):
    from tumult.synthesis import submission, workflows
    request = _request(Atmosphere1D('one.mod', (5750, 4.5, 0, 1)), tmp_path)

    def execute(resolved, cancel_token=None):
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), ({'window': None},), resolved.request,
                               resolved.settings)

    monkeypatch.setattr(workflows, 'execute_request', execute)
    handle = submission.submit(request, root=tmp_path / 'jobs')
    result = handle.result(timeout=5)
    reopened = submission.open_workflow(handle.folder)
    assert reopened.status == 'finished'
    assert np.allclose(reopened.result().spectrum.values[1], result.spectrum.values[1])
    assert reopened.progress['phase'] == 'finished'


def test_submission_records_scheduled_runs_before_completion(tmp_path, monkeypatch,
                                                            workflow_daemon):
    """The workflow record names its native runs as the executor schedules them."""
    import threading
    from tumult.synthesis import runtime, submission, workflows
    request = _request(Atmosphere1D('one.mod', (5750, 4.5, 0, 1)), tmp_path)
    scheduled = tmp_path / 'out' / 'run_a'
    release = threading.Event()
    seen = {}

    def execute(resolved, cancel_token=None):
        runtime.note_run(str(scheduled))            # what _run_batch does per job
        seen['during'] = submission.open_workflow(handle.folder).runs
        release.wait(5)
        spec = Spectrum(np.array([5000., 5001.]), np.array([1., .9]))
        return SynthesisResult((spec,), ({'window': None},), resolved.request,
                               resolved.settings, runs=(str(scheduled),))

    monkeypatch.setattr(workflows, 'execute_request', execute)
    handle = submission.submit(request, root=tmp_path / 'jobs')
    release.set()
    handle.result(timeout=5)
    assert seen['during'] == (str(scheduled),)
    assert submission.open_workflow(handle.folder).runs == (str(scheduled),)


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
        "&composition_params\n abundances='Mg=0.2'\n/\n")
    return Departures.load(folder)


def test_check_request_refuses_unsupported_combinations_before_a_record_exists(tmp_path):
    from tumult.synthesis import check_request, submission
    prefs = {'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out', 'R': 50000}
    field = _saved_departures(tmp_path)
    one = Atmosphere1D('one.mod', (5750, 4.5, 0, 1))
    window = [MaskWindow((5000, 5001), (5000, 5001))]
    supported = NativeRequest(prefs, run_vals={'mode': 'NLTE', 'atoms': ['mg.atom']},
                             selection=one, masks=window, departures=field)
    assert check_request(supported) is supported
    rejected = [
        (NativeRequest(prefs, run_vals={'mode': 'LTE'}, selection=one, masks=window,
                      departures=field), 'NLTE mode'),
        (NativeRequest(prefs, run_vals={'mode': 'NLTE', 'atoms': ['mg.atom']},
                      selection=GridTarget({'folder': ''}, (5750, 4.5, 0, 1)), masks=window,
                      departures=field), 'DepartureGrid'),
        (NativeRequest(prefs, run_vals={'mode': 'NLTE', 'atoms': []}, selection=one,
                      masks=window), 'model atoms'),
    ]
    for request, cause in rejected:
        with pytest.raises(ValueError, match=cause):
            check_request(request)
        with pytest.raises(ValueError, match=cause):
            request.resolve()
        with pytest.raises(ValueError, match=cause):
            submission.submit(request, root=tmp_path / 'jobs')
    assert not (tmp_path / 'jobs').exists()      # nothing was recorded for a refused request


def test_fit_request_refuses_what_fitting_does_not_support(tmp_path):
    from tumult.synthesis import FitRequest
    prefs = {'dispatch_cwd': str(tmp_path), 'R': 50000}
    masks = [MaskWindow((5000.2, 5000.8), (5000, 5001)),
             MaskWindow((5002.2, 5002.6), (5002, 5003))]
    grid = NativeRequest(prefs, run_vals={'mode': 'LTE', 'abundances': ''},
                        selection=GridTarget({'folder': 'grid'}, (5750, 4.5, 0, 1)),
                        masks=masks)
    fit = FitRequest(grid, fit_stellar=['teff'], abund_fits={'Mg': 0.1},
                     fit_post={'vshift': 0})
    assert fit.synth_span == (5000., 5003.)               # the span every trial synthesises
    assert fit.fit_windows == ((5000.2, 5000.8), (5002.2, 5002.6))   # compared inside
    assert FitRequest.from_dict(fit.to_dict()).to_dict() == fit.to_dict()
    assert fit.axes == ['teff', 'Mg']

    one = Atmosphere1D('one.mod', (5750, 4.5, 0, 1))
    with pytest.raises(ValueError, match='1D atmospheres and grids'):
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'LTE'}, masks=masks,
                                selection=Snapshot3D('cube.dat', (5750, 4.5, 0), (4, 4, 8))),
                   fit_post={'vshift': 0})
    with pytest.raises(ValueError, match='saved departures'):
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'NLTE', 'atoms': ['mg.atom']},
                                selection=one, masks=masks,
                                departures=_saved_departures(tmp_path)),
                   abund_fits={'Mg': 0.0})
    with pytest.raises(ValueError, match='single atmosphere file'):
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'LTE'}, selection=one, masks=masks),
                   fit_stellar=['teff'])
    with pytest.raises(ValueError, match='at least one parameter'):
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'LTE'}, selection=one, masks=masks))
    with pytest.raises(ValueError, match='window'):
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'LTE'}, selection=one),
                   abund_fits={'Mg': 0.0})
    with pytest.raises(ValueError, match='model atoms'):    # check_request applies too
        FitRequest(NativeRequest(prefs, run_vals={'mode': 'NLTE'}, selection=one, masks=masks),
                   abund_fits={'Mg': 0.0})


def test_fit_executes_each_trial_point_as_the_request_over_the_span(tmp_path, monkeypatch):
    """Every trial point is the fit's own request at that point, over one span window, so
    preferences, expert settings, sampling and padding reach dispatch unchanged; the
    observation is compared only inside the fit windows."""
    from tumult.synthesis import FitRequest, workflows
    prefs = {'dispatch_cwd': str(tmp_path), 'R': 50000, 'abund_file': None,
             'abund_step': 0.5, 'fit_max_rounds': 3}
    masks = [MaskWindow((5000.2, 5000.8), (5000, 5001)),
             MaskWindow((5002.2, 5002.6), (5002, 5003))]
    synthesis = NativeRequest(prefs, run_vals={'mode': 'LTE', 'abundances': 'Ca=0.1',
                                              'relative_abundances': True},
                             selection=GridTarget({'folder': 'grid', 'format': 'Marcs'},
                                                  (5750, 4.5, 0, 1)),
                             masks=masks, expert_settings={'spectrum_params': {'daa': 0.01}})
    request = FitRequest(synthesis, fit_stellar=['teff'], abund_fits={'Mg': 0.0},
                         fit_post={'vshift': 0.0})
    executed = []

    def execute(varied, cancel_token=None):
        executed.append(varied)
        xx = np.linspace(5000., 5003., 31)
        return SynthesisResult((Spectrum(xx, np.ones_like(xx)),),
                               ({'window': [5000., 5003.]},), varied, {})

    nodes = {(5500., 4.5, 0., 1.): 'a', (6000., 4.5, 0., 1.): 'b'}
    monkeypatch.setattr(workflows, 'execute_request', execute)
    monkeypatch.setattr(workflows, 'discover_grid', lambda grid, root=None: nodes)
    seen = {}

    def adaptive(target, synth_point_fn, axes, guesses, **kwargs):
        seen.update(target=target, axes=axes, guesses=guesses, **kwargs)
        for teff in (5500., 6000.):
            synth_point_fn({'teff': teff, 'Mg': 0.25})
        best = {'teff': 5600., 'Mg': 0.25, 'vshift': 1.5}
        return best, synth_point_fn(best), None

    monkeypatch.setattr(workflows, 'adaptive_fit', adaptive)
    xx = np.arange(4999., 5004., 0.1)
    observation = Spectrum(xx, np.ones_like(xx))
    best, raw, post = workflows.fit(observation, request)

    assert best['teff'] == 5600. and post == {'vshift': 1.5}
    assert seen['axes'] == ['teff', 'Mg'] and seen['post_fit'] == ['vshift']
    assert seen['axis_nodes'] == {'teff': [5500., 6000.]} and seen['spacings'] == {'Mg': 0.5}
    assert seen['max_rounds'] == 3
    inside = seen['target'].values[0]
    assert inside.min() >= 5000.2 and inside.max() <= 5002.6
    assert not ((inside > 5000.8) & (inside < 5002.2)).any()   # the gap is dropped
    # the trial points (and the raw recompute of the best one): the request at each point
    # over the one span window
    trial = executed[0]
    assert [v.selection.stellar[0] for v in executed] == [5500., 6000., 5600., 5600.]
    assert trial.run_vals['teff'] == 5500. and trial.masks == (MaskWindow((5000., 5003.),
                                                                          (5000., 5003.)),)
    assert trial.run_vals['abundances'] == 'Ca=0.1 Mg=0.2500'
    assert trial.prefs == synthesis.prefs
    assert trial.expert_settings == synthesis.expert_settings
    assert type(trial.selection) is GridTarget
    assert trial.selection.grid == synthesis.selection.grid
    assert raw.values[0][0] == 5000. and raw.values[0][-1] == 5003.


@pytest.mark.parametrize('with_linelist', [False, True])
def test_snapshot_correction_uses_atom_grid_without_unbounded_resolving_power(
        tmp_path, with_linelist):
    from tumult.synthesis import snapshot

    atom = tmp_path / 'atom.ba'
    line_grid = tmp_path / 'line_grid.txt'
    request = NativeRequest(
        {'dispatch_cwd': str(tmp_path), 'R': 1e6, 'linelist_crop_aa': 0},
        run_vals={'mode': 'NLTE', 'atoms': [str(atom)], 'atom_lines': True,
                  'only_atom': not with_linelist, 'line_grid': str(line_grid)},
        selection=Snapshot3D('cube.atmos', (4548, 3.964, 0), (20, 20, 64)))
    resolved = request.resolve()
    field = snapshot.Departures(
        folder=tmp_path / 'departures', name=Path('synth/dep'),
        atmosphere='cube.atmos', resolution=(20, 20, 64),
        settings=snapshot.departure_settings(resolved.settings, 'cube.atmos'))
    formal = snapshot.spectrum_settings(field, resolved.settings)
    assert not formal.get('spectrum_params', {}).get('R')
    assert formal['atom_params']['line_grid_file'] == str(line_grid)
    assert formal['atmos_params']['read_departures']
    assert formal['rt_iteration_params']['maxiter'] == 0

    # A user-selected wavelength window still needs the requested sampling.
    windowed = request.with_masks([MaskWindow((4553, 4555), (4552, 4556))]).resolve()
    assert windowed.settings['spectrum_params']['R'] == 1e6
