"""Pin the broadening kernels (Gray 1992 rotation, Gaussian, radial-tangential macro)
and the full Spectrum.crop convolution path the GUI depends on."""
import astropy.units as u
import numpy as np
import pytest

from cases import convolution_cases, gauss_line
from conftest import assert_matches_golden, load_golden

kms = u.km / u.s


def test_matches_golden():
    assert_matches_golden(convolution_cases(), load_golden("convolutions"))


def test_kernels_are_normalised():
    from tumult.convolutions import combined_kernel

    px, py = combined_kernel(rd=2e-6, vdop=3.0 * kms, vrot=2.0 * kms, zeta_rt=1.5 * kms)
    area = np.trapz(np.asarray(py, float), x=np.asarray(px, float))
    assert area == pytest.approx(1.0, rel=1e-6)


def test_convolution_conserves_equivalent_width():
    """Broadening redistributes flux; it must conserve the line's EW (physics invariant)."""
    from tumult.spectrum import Spectrum

    x, y = gauss_line(n=4001, span=4.0)
    sp = Spectrum(x, y)
    before = sp.get_weq(x0=5000.0 * u.AA, profile="trapz", return_profile=True).mAA
    conv = sp.crop(vdop=2.0, vrot=1.5, zeta_rt=1.0)
    after = conv.get_weq(x0=5000.0 * u.AA, profile="trapz", return_profile=True).mAA
    assert after == pytest.approx(before, rel=1e-3)


def test_vshift_is_a_pure_doppler_shift():
    """crop(vshift=v) keeps the original wavelength grid and shifts the *flux*:
    y_out(x) = y(x / (1+beta)), implemented as interp(x, x*(1+beta), y)."""
    from tumult.spectrum import Spectrum, get_wavelength_shift

    x, y = gauss_line()
    sp = Spectrum(x, y)
    shifted = sp.crop(vshift=10.0)
    factor = float(get_wavelength_shift(10.0))
    np.testing.assert_allclose(np.asarray(shifted.xdata.value, float), x, rtol=0)
    np.testing.assert_allclose(
        np.asarray(shifted.ydata, float),
        np.interp(x, x * factor, y), rtol=1e-12,
    )


def test_custom_kernel_call_is_broken():
    """DOCUMENTS A LIVE BUG (found in review, to fix in Phase 1): combined_kernel calls
    custom_kernel(custom, rd, px) but the function signature is custom_kernel(custom, px).
    Any custom-kernel convolution currently raises TypeError. When Phase 1 fixes it,
    replace this test with a golden for the fixed behaviour."""
    from tumult.convolutions import combined_kernel

    xc = np.linspace(-5, 5, 11) * kms
    yc = np.exp(-0.5 * (np.linspace(-5, 5, 11) / 2.0) ** 2)
    with pytest.raises(TypeError):
        combined_kernel(rd=2e-6, custom=(xc, yc))
