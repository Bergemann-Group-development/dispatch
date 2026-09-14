"""Interaction contracts that change scientific state or preserve user work."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'TumultStudio2'))

import store
from controls import Field, Form
from draft import CalculationDraft
from toolmanager import ToolManager
from tumult import Spectrum
from tumult.linelist import load_master
from tumult.synthesis import Atmosphere1D, NativeRequest, MaskWindow
from tumult.synthesis import workflows


def test_cancel_draft_does_not_persist_and_restores_fields(tmp_path, monkeypatch):
    monkeypatch.setattr(store, 'CONFIG_PATH', tmp_path / 'prefs.json')
    form = Form('measurement-test', [Field('snr', 'SNR', value=100)])
    draft = CalculationDraft(forms=[('measurement', form)])
    draft.begin()
    form.set_values(snr=25)
    assert store.get_section('measurement') == {}
    draft.discard()
    assert form.value('snr') == 100
    assert not store.CONFIG_PATH.exists()
    draft.begin()
    form.set_values(snr=50)
    draft.commit()
    draft.discard()  # successful submission's close event must not roll back
    assert store.get_section('measurement')['snr'] == form.value('snr') == 50


def test_unavailable_tool_preserves_current_view():
    current = SimpleNamespace(active=True, deactivate=lambda: pytest.fail('lost current view'))
    messages = []
    unavailable = SimpleNamespace(active=False, canvas=SimpleNamespace(selected=lambda: None),
                                  notify=messages.append)
    tools = ToolManager()
    tools.register('linemask', current)
    tools.register('normalise', unavailable)
    tools.activate('normalise')
    assert tools.active_name() == 'linemask'
    assert messages


def test_line_list_folder_loads_atomic_files_and_hydrogen(tmp_path):
    """The configured line list may be the synthesis folder: navigation reads its atomic
    files and HOP's dedicated hydrogen file, and leaves the molecular lists to HOP."""
    for name, species, lam in [('fe', 'Fe I', 5000), ('ch', 'CH', 5001)]:
        (tmp_path / f'{name}.txt').write_text(
            f"'{species}'\n'0' 1\n{lam} 1 -1 1 1 1\n")
    (tmp_path / 'Hlinedata').write_text(
        "'H I  lambda(AA) nblo nbup elo(eV) eup(eV) loggf ident raddamp'\n'1' 1\n"
        "6562.797 2 3 10.199 12.088 0.710 'Ha' 5.70E+08\n")
    frame = load_master(tmp_path)
    assert set(frame.ion) == {'Fe I', 'H I'}
    assert frame[frame.ion == 'H I'].iloc[0][['lam0', 'evi', 'loggf']].tolist() == [6562.797, 10.199, 0.71]
    assert frame.attrs['strength_is_proxy']


def test_atom_only_request_skips_background_preparation(tmp_path, monkeypatch):
    request = NativeRequest({'dispatch_cwd': str(tmp_path), 'linelist': 'background',
                            'linelist_crop_aa': 5},
        run_vals={'mode': 'NLTE', 'atoms': ['Fe.atom'], 'only_atom': True},
        selection=Atmosphere1D('science.mod', (5750, 4.5, 0, 1)),
        masks=[MaskWindow((5000, 5001), (5000, 5001))])
    resolved = request.resolve()
    assert resolved.settings['linelist_params'] == {'linelist_file': [''], 'linelist_folder': ''}
    assert workflows._prepare_artifact(resolved) is None
    assert NativeRequest.from_json(request.to_json()).run_vals['only_atom'] is True


def test_async_add_preserves_selection_and_wavelength_view():
    from canvas import SpectrumCanvas
    canvas = SpectrumCanvas()
    first = canvas.add_object(Spectrum(np.linspace(5000, 5010, 100), np.ones(100)))
    canvas.main.x_range.start, canvas.main.x_range.end = 5002, 5003
    canvas.add_object(Spectrum(np.linspace(6000, 6010, 100), np.ones(100)),
                      source=('file', 'second.txt'), select=False)
    assert canvas.selected() is first
    assert (canvas.main.x_range.start, canvas.main.x_range.end) == (5002, 5003)


def test_unusable_line_list_never_reaches_hop(tmp_path):
    """The configured list is a user preference now, so an empty, missing or foreign path
    must be answered by the proxy source; HOP aborts the process on such input."""
    from tumult.linestrength import LineStrengths, usable_linelist
    (tmp_path / 'notes.txt').write_text('no line list here\n')
    assert not usable_linelist('') and not usable_linelist(tmp_path / 'missing')
    assert not usable_linelist(tmp_path)
    good = tmp_path / 'fe.txt'
    good.write_text("'Fe I'\n'0' 1\n5000 1 -1 1 1 1\n")
    assert usable_linelist(good) and usable_linelist(tmp_path)
    strengths = LineStrengths(master_path=tmp_path / 'missing', linelist=tmp_path / 'missing')
    obj = SimpleNamespace(meta={}, source=None)
    assert not strengths.pending(obj, 5000., 5001.)
    assert strengths.warm(obj, 5000., 5001.).lines is None
