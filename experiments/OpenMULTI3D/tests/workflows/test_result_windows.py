"""SynthesisResult.crop and plot_line_labels choose the windows from the wavelength selection.

The sources are stand-ins for the spectrum runs of two windows, 5000-5010 and 5020-5030 Å,
whose spectra are their wavelengths, so a crop shows which window and range it came from.
"""
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tumult.spectrum import Spectrum
from tumult.synthesis import SynthesisResult

WINDOWS = ((5000.0, 5010.0), (5020.0, 5030.0))


class Source:
    """A run over one window: crops of y = x on a 0.1 Å grid, and labels at its middle."""
    def __init__(self, bounds):
        self.xx = np.linspace(bounds[0] - 1, bounds[1] + 1, int((bounds[1] - bounds[0] + 2) * 10) + 1)
        self.run = self
        self.middle = sum(bounds) / 2

    def crop(self, LTE=False, **kwargs):
        return Spectrum(self.xx, self.xx).crop(**kwargs)

    def plot_line_labels(self, xmin, xmax, **kwargs):
        lam0 = [lam for lam in (self.middle,) if xmin < lam < xmax]
        return pd.DataFrame({"lam0": lam0})


@pytest.fixture
def result():
    spectra = tuple(Source(bounds).crop(xmin=bounds[0], xmax=bounds[1]) for bounds in WINDOWS)
    return SynthesisResult(spectra, WINDOWS, SimpleNamespace(is_nlte=False), {}, {},
                           sources=tuple(Source(bounds) for bounds in WINDOWS))


def test_without_a_selection_the_crop_covers_every_window(result):
    spec = result.crop()
    assert spec.n_segments == 2
    for part, window in zip(spec.split_segments(), WINDOWS):
        assert part.xdata.value[[0, -1]].tolist() == pytest.approx(list(window))


def test_a_selection_inside_one_window_is_an_ordinary_spectrum(result):
    for spec in (result.crop(x0=5025 * u.AA, dx=1), result.crop(x0=5025, width=2),
                 result.crop(window=(5024, 5026)), result.crop(xmin=5024, xmax=5026)):
        assert spec.segments is None
        assert spec.xdata.value[[0, -1]].tolist() == pytest.approx([5024, 5026])


def test_an_open_selection_ends_at_the_window_bound(result):
    spec = result.crop(xmin=5025)
    assert spec.xdata.value[[0, -1]].tolist() == pytest.approx([5025, 5030])


def test_a_selection_across_windows_keeps_each_part(result):
    spec = result.crop(window=(5005, 5025))
    blue, red = spec.split_segments()
    assert blue.xdata.value[[0, -1]].tolist() == pytest.approx([5005, 5010])
    assert red.xdata.value[[0, -1]].tolist() == pytest.approx([5020, 5025])
    newx = np.array([5001.0, 5002.0, 5021.0])
    assert result.crop(newx=newx).xdata.value.tolist() == newx.tolist()


def test_a_selection_outside_every_window_is_refused(result):
    with pytest.raises(ValueError, match="overlaps no window"):
        result.crop(window=(5012, 5018))
    with pytest.raises(ValueError, match="inside one window"):
        result.crop(slc=slice(0, 5))


def test_labels_come_from_every_window_in_range(result):
    assert result.plot_line_labels(xmin=5000, xmax=5030)["lam0"].tolist() == [5005, 5025]
    assert result.plot_line_labels(xmin=5015, xmax=5030)["lam0"].tolist() == [5025]
