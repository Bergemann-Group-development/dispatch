"""Segmented output grids: reading ``out_segments.txt`` and respecting the gaps.

A multi-window run concatenates disjoint windows into one x array, so the gaps are not
recoverable from the wavelengths alone and any neighbourhood operation across one mixes
values that are not neighbours. These tests pin the two things that keeps honest: the
segment table round-trips to the windows that were asked for, and an operation needing a
neighbourhood refuses a spectrum that still spans several of them.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from tumult.flux import read_segments
from tumult.io_binary import write_patch_save
from tumult.run import Run
from tumult.spectrum import Spectrum
from tumult.synthesis import compute_windows, normalise_windows
from tumult.synthesis.api import window_file_text

HEADER = "# segment  first  last  blue_vac_AA  red_vac_AA\n"


def write_segments(folder, rows):
    save = folder / 'save'
    save.mkdir(parents=True, exist_ok=True)
    (save / 'out_segments.txt').write_text(HEADER + "".join(rows))
    return read_segments(folder)


def test_no_file_means_contiguous(tmp_path):
    """A single-window run writes no table and needs none."""
    (tmp_path / 'save').mkdir()
    assert read_segments(tmp_path) is None


def test_segments_round_trip_to_air_windows(tmp_path):
    """Indices survive as given; bounds come back in air, matching Run.lam."""
    segs = write_segments(tmp_path, [
        "        1      1    195    5151.43473    5171.44006\n",
        "        2    196    389    5181.44272    5201.44805\n",
    ])
    assert [s['first'] for s in segs] == [1, 196]
    assert [s['last'] for s in segs] == [195, 389]
    # the vacuum bounds were built from air windows 5150-5170 and 5180-5200
    assert segs[0]['blue'] == pytest.approx(5150.0, abs=1e-3)
    assert segs[0]['red'] == pytest.approx(5170.0, abs=1e-3)
    assert segs[1]['blue'] == pytest.approx(5180.0, abs=1e-3)
    assert segs[1]['red'] == pytest.approx(5200.0, abs=1e-3)



def test_windows_split_into_several_patches_stitch_to_one_spectrum(tmp_path):
    """With fewer windows than threads each window has several output patches, saved as
    flux_008001..; all of them are the output spectrum, whatever the internal n_windows."""
    segs = write_segments(tmp_path, [
        "        1      1      8    5151.43473    5171.44006\n",
        "        2      9     16    5181.44272    5201.44805\n",
    ])
    flux = np.arange(16, dtype=np.float32)
    parts = np.split(flux, 8)   # four patches per window
    write_patch_save(tmp_path / 'save', 'flux',
                     {(8, i + 1): p.reshape(-1, 1) for i, p in enumerate(parts)})
    run = SimpleNamespace(nml={'run_params': {'n_windows': 4}}, segments=segs, sfolder=tmp_path / 'save')
    out = Run.read_patch_save(run, 'flux', lazy=False)
    assert len(out) == 1
    np.testing.assert_array_equal(np.asarray(out[0]).ravel(), flux)

def two_window_spectrum():
    """Two disjoint windows concatenated, as a multi-window run hands them over."""
    x = np.concatenate([np.linspace(5150, 5170, 21), np.linspace(5180, 5200, 21)])
    y = np.concatenate([np.full(21, 1.0), np.full(21, 2.0)])
    segments = [{'first': 1, 'last': 21, 'blue': 5150.0, 'red': 5170.0},
                {'first': 22, 'last': 42, 'blue': 5180.0, 'red': 5200.0}]
    return Spectrum(x, y, segments=segments)


def test_split_segments_slices_by_index():
    """Endpoints sit exactly on their bounds, so the split must not go via floats."""
    parts = two_window_spectrum().split_segments()
    assert [len(p.xdata) for p in parts] == [21, 21]
    assert [p.n_segments for p in parts] == [1, 1]
    assert np.asarray(parts[0].xdata.value)[-1] == pytest.approx(5170.0)
    assert np.asarray(parts[1].xdata.value)[0] == pytest.approx(5180.0)


def test_contiguous_spectrum_splits_to_itself():
    """Callers can loop over split_segments() without special-casing one window."""
    sp = Spectrum(np.linspace(5150, 5170, 21), np.ones(21))
    assert sp.n_segments == 1
    assert sp.split_segments() == [sp]


def test_convolution_across_a_gap_is_refused():
    sp = two_window_spectrum()
    with pytest.raises(ValueError, match="contiguous"):
        sp.convolve(res=50000)


def test_measurement_across_a_gap_is_refused():
    sp = two_window_spectrum()
    with pytest.raises(ValueError, match="contiguous"):
        sp.measure()


def test_crop_into_one_window_is_an_ordinary_spectrum():
    """This is how a caller gets from a multi-window run to something measurable."""
    sub = two_window_spectrum().crop(window=(5152, 5168))
    assert sub.n_segments == 1
    assert sub.convolve(res=50000) is not None


# --- the multi-window synthesis API -------------------------------------------------


def test_windows_are_sorted_and_padded():
    """Order of request must not matter, and each window gets convolution margin."""
    assert normalise_windows([(5180, 5200), (5150, 5170)], pad=1.0) == [
        (5149.0, 5171.0), (5179.0, 5201.0)]


def test_touching_windows_merge():
    """Overlapping requests would otherwise make segments sharing grid points."""
    assert normalise_windows([(5150, 5170), (5169, 5185)], pad=1.0) == [(5149.0, 5186.0)]


def test_reversed_bounds_are_accepted():
    assert normalise_windows([(5170, 5150)], pad=1.0) == [(5149.0, 5171.0)]


def test_default_padding_scales_with_width():
    """Matches _apply_window's rule, so one window here covers what compute_spectrum does."""
    (blue, red), = normalise_windows([(5000, 5100)])
    assert (blue, red) == (4995.0, 5105.0)


@pytest.mark.parametrize("bad", [[], [(0.0, 10.0)], [(np.nan, 5000.0)]])
def test_invalid_windows_are_rejected(bad):
    with pytest.raises(ValueError):
        normalise_windows(bad)


def test_window_file_text_is_deterministic():
    """The cache keys on this text, so equal requests must produce equal bytes."""
    a = window_file_text(normalise_windows([(5150, 5170), (5180, 5200)], pad=1.0))
    b = window_file_text(normalise_windows([(5180, 5200), (5150, 5170)], pad=1.0))
    assert a == b
    assert a.splitlines()[1:] == ["5149.0000  5171.0000", "5179.0000  5201.0000"]


def test_compute_windows_requires_R():
    """A window file fixes per-window sampling; daa/nnu would leave the grid ambiguous."""
    with pytest.raises(ValueError, match="positive R"):
        compute_windows('unused', (5777, 4.4, 0.0), {}, windows=[(5150, 5170)], R=None)
