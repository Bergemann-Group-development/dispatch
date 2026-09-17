"""The extracted line workflows run without the GUI, with the synthesis stubbed out."""
import numpy as np
import pytest

from tumult import abundance as ab
from tumult.abundance import LineResult
from tumult.linemask import LineMask
from tumult.spectrum import Spectrum


def _spectrum(depth):
    """A Gaussian line whose depth (hence EW) grows with the abundance offset."""
    x = np.arange(4999.0, 5001.0, 0.01)
    y = 1.0 - depth * np.exp(-0.5 * ((x - 5000.0) / 0.08) ** 2)
    return Spectrum(x, y)


def _mask():
    return LineMask(lam0=5000.0, lmin=4999.4, lmax=5000.6, smin=4999.2, smax=5000.8)


def test_fit_windows_restrict_the_target_samples():
    """A line-mask fit compares only inside the fit windows: the samples between them are
    dropped from the target, while the synthesised span stays continuous."""
    from tumult.synthesis.operations import _inside_windows

    x = np.arange(5000., 5010., 0.5)
    spec = Spectrum(x, np.ones_like(x))
    sub = _inside_windows(spec, [(5001.0, 5002.0), (5008.0, 5009.0)])
    assert list(sub.values[0]) == [5001.0, 5001.5, 5002.0, 5008.0, 5008.5, 5009.0]
    assert sub.xunit == spec.xunit and sub.yunit == spec.yunit


def test_profile_fit_uses_fixed_model_continuum():
    """A fixed local continuum multiplies the model while fitting; the overlay stays
    normalised and the target is left untouched."""
    nodes = ((4999.0, 1.02), (5001.0, 1.06))
    line = LineResult(lam0=5000.0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=90.0, ew_err=3.0,
                      fwhm=0.19, mask=_mask(), renorm=True, continuum=nodes)
    target = _spectrum(0.35 * 10 ** (0.5 * 0.23)).crop(continuum=list(nodes))
    original = target.ydata.copy()

    def synthesize(offsets):
        return [[_spectrum(0.35 * 10 ** (0.5 * a))] for a in offsets]

    ab.fit_profiles([line], synthesize, target, step=0.25, limit=1.0, max_rounds=8,
                    profile="gauss")
    assert line.a == pytest.approx(0.23, abs=0.02)
    assert line.fit_spec.ydata.value.max() == pytest.approx(1.0, abs=1e-5)
    np.testing.assert_array_equal(target.ydata, original)
