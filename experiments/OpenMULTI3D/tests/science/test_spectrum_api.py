"""Spectrum.measure, Spectrum.convolve and crop(window=...) against the older spellings.

The older spellings (get_weq flags, crop's four window styles, crop(**conv)) remain in
use and are pinned by goldens. These tests show that both forms give identical results.
"""
import astropy.units as u
import numpy as np
import pytest

from tests.support.spectra import LINE_LAM0 as LINE_X0, gaussian_line


def _spectrum():
    from tumult.spectrum import Spectrum

    x, y = gaussian_line()
    return Spectrum(x, y), x, y


def test_measure_equals_get_weq_across_profiles():
    from tumult.profile import profile_names

    sp, x, y = _spectrum()
    for profile in profile_names():
        via_measure = sp.measure(profile=profile, x0=LINE_X0)  # plain x0: coerced
        via_get_weq = sp.get_weq(profile=profile, x0=LINE_X0 * u.AA, return_profile=True)
        assert via_measure.maa == pytest.approx(via_get_weq.maa, rel=1e-12), profile
        assert via_measure.rew == pytest.approx(via_get_weq.rew, rel=1e-12), profile


def test_measure_window_equals_manual_crop():
    sp, x, y = _spectrum()
    window = (4999.0, 5001.0)
    direct = sp.measure(window, profile="trapz", x0=LINE_X0)
    manual = sp.crop(xmin=window[0], xmax=window[1]).measure(profile="trapz", x0=LINE_X0)
    assert direct.maa == pytest.approx(manual.maa, rel=1e-12)


def test_crop_window_equals_xmin_xmax():
    sp, x, y = _spectrum()
    a = sp.crop(window=(4999.0, 5001.0))
    b = sp.crop(xmin=4999.0, xmax=5001.0)
    np.testing.assert_array_equal(np.asarray(a.xdata.value), np.asarray(b.xdata.value))
    np.testing.assert_array_equal(np.asarray(a.ydata), np.asarray(b.ydata))


def test_crop_window_mutually_exclusive():
    sp, x, y = _spectrum()
    with pytest.raises(ValueError, match="mutually exclusive"):
        sp.crop(window=(4999.0, 5001.0), xmin=4999.0)
    with pytest.raises(ValueError, match="mutually exclusive"):
        sp.crop(window=(4999.0, 5001.0), x0=5000.0, dx=1.0)


def test_convolve_equals_crop_conv():
    sp, x, y = _spectrum()
    a = sp.convolve(vdop=3.0, vrot=2.0, zeta_rt=1.5, vshift=0.5)
    b = sp.crop(vdop=3.0, vrot=2.0, zeta_rt=1.5, vshift=0.5)
    np.testing.assert_array_equal(np.asarray(a.xdata.value), np.asarray(b.xdata.value))
    np.testing.assert_array_equal(np.asarray(a.ydata, float), np.asarray(b.ydata, float))
    assert sp.convolve() is sp  # no-op form


def test_fit_line_profile_method_fixed():
    """The Spectrum.fit_line_profile method used to unpack 4 values from a 3-tuple and
    could never run; it now returns the fitted model as a Spectrum."""
    from tumult.spectrum import Spectrum

    sp, x, y = _spectrum()
    fitted = sp.fit_line_profile(profile="gauss", x0=LINE_X0)
    assert isinstance(fitted, Spectrum)
    # the fitted gaussian should reproduce the noiseless input almost exactly
    np.testing.assert_allclose(np.asarray(fitted.ydata, float), y, atol=1e-8)
