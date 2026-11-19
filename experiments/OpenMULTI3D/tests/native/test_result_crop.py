"""A synthesis result crops its NLTE and LTE spectra again from the runs it came from.

An NLTE spectrum run also computes the LTE spectrum, so :meth:`SynthesisResult.crop` gives
both without another run, for an individual atmosphere and for a saved and reopened
result. The case is Ba II 4554 Å in the MARCS solar model with the six-level toy atom of
``examples/sun_ba_1d.nml``, without a line list (a few seconds). The result also labels
the lines of its run.

Needs the built ``dispatch.x``.
"""
from pathlib import Path

import numpy as np
import pytest

import tumult
from tumult import departures, synthesis
from tumult.synthesis import Star, SynthConfig, SynthesisResult, SynthRequest

ROOT = Path(__file__).resolve().parents[2]
WINDOW = (4552.0, 4556.0)

pytestmark = pytest.mark.skipif(not (ROOT / 'dispatch.x').is_file(), reason='needs dispatch.x')


@pytest.fixture(scope='module')
def result(tmp_path_factory):
    root = tmp_path_factory.mktemp('crop')
    for name in ('dispatch.x', 'input'):
        (root / name).symlink_to(ROOT / name)
    config = SynthConfig(root=root, cache='synth', abundances='input/abund/abund_magg22.txt',
                         threads=2)
    request = SynthRequest(config=config, atmosphere='input/atmos/sun.mod',
                           star=Star(5777, 4.44, 0.0, 1.0), windows=WINDOW)
    nlte = departures.request(request, atom='input/atoms/atom.ba06')
    return synthesis.synthesize(request.update(departures=nlte))


def _same(a, b):
    np.testing.assert_array_equal(a.xdata.value, b.xdata.value)
    np.testing.assert_array_equal(np.asarray(a.ydata), np.asarray(b.ydata))


def test_crop_gives_the_result_spectrum(result):
    _same(result.crop(), result.spectrum)


def test_lte_crop_is_the_lte_spectrum_of_the_run(result):
    run = tumult.read(result.runs[-1])
    lte = result.crop(LTE=True)
    _same(lte, run.crop(xmin=WINDOW[0], xmax=WINDOW[1], LTE=True))
    nlte_core, lte_core = result.crop().ydata.min(), lte.ydata.min()
    assert abs(nlte_core - lte_core) > 0.01


def test_result_labels_the_lines_of_its_run(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    result.crop().plot()
    labels = result.plot_line_labels(show_atom=True)
    plt.close('all')
    assert any(abs(lam - 4554.03) < 0.05 for lam in labels['lam0'])


def test_a_reopened_result_crops_the_same(result, tmp_path):
    reopened = SynthesisResult.load(result.save(tmp_path / 'saved'))
    for lte in (False, True):
        _same(reopened.crop(LTE=lte), result.crop(LTE=lte))
