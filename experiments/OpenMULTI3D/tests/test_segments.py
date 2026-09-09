"""Segmented output grids: reading ``out_segments.txt`` and respecting the gaps.

A multi-window run concatenates disjoint windows into one x array, so the gaps are not
recoverable from the wavelengths alone and any neighbourhood operation across one mixes
values that are not neighbours. These tests pin the two things that keeps honest: the
segment table round-trips to the windows that were asked for, and an operation needing a
neighbourhood refuses a spectrum that still spans several of them.
"""
import numpy as np
import pytest

from tumult.flux import read_segments
from tumult.spectrum import Spectrum

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
