"""Public request, preparation, departure and durable-submission contracts."""
import time
from pathlib import Path

import numpy as np

from tumult import Spectrum
from tumult.synthesis import (Atmosphere1D, GridTarget, MaskWindow, ResolvedSynthRequest, Snapshot3D,
                              SynthesisResult, SynthRequest)


def _request(selection, tmp_path, masks=()):
    return SynthRequest({'dispatch_cwd': str(tmp_path), 'synth_output_dir': 'out', 'R': 50000},
                        run_vals={'mode': 'LTE'}, selection=selection, masks=masks)


def test_all_atmosphere_selections_roundtrip_and_capture_masks(tmp_path):
    masks = [{'fit': [5000, 5001], 'synth': [4999.5, 5001.5], 'include': True}]
    selections = [GridTarget({'folder': 'grid'}, (5750, 4.5, 0, 1)),
                  Atmosphere1D('one.mod', (5750, 4.5, 0, 1)),
                  Snapshot3D('cube.dat', (5750, 4.5, 0), (8, 8, 16))]
    for selection in selections:
        request = _request(selection, tmp_path, masks)
        masks[0]['synth'][0] = 1
        reopened = SynthRequest.from_json(request.to_json())
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

    def crop(path, *_args, **_kwargs):
        calls.append(Path(path))
        Path(path).write_text('prepared')

    monkeypatch.setattr(linestrength, 'crop_linelist', crop)
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


def test_departure_field_can_be_opened_from_native_artifacts(tmp_path):
    from tumult.synthesis import DepartureField
    folder = tmp_path / 'dep'
    (folder / 'save').mkdir(parents=True)
    (folder / 'save' / 'pop_Mg_meta.txt').write_text('')
    (folder / 'params.nml').write_text(
        "&atmos_params\n atmos_file='cube.dat' nx=4 ny=5 nz=6\n/\n"
        "&composition_params\n abundances='Mg=0.2'\n/\n")
    field = DepartureField.load(folder)
    assert field.saved and field.resolution == (4, 5, 6)
    assert field.species == ('Mg',) and field.composition['abundances'] == 'Mg=0.2'


def test_durable_submission_reopens_saved_result(tmp_path, monkeypatch):
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
