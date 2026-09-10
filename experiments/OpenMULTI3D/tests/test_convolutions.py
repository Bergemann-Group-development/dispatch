"""Pin the broadening kernels (Gray 1992 rotation, Gaussian, radial-tangential macro)
and the full Spectrum.crop convolution path the GUI depends on."""
import astropy.units as u
import numpy as np
import pytest

from cases import convolution_cases, gauss_line
from conftest import assert_matches_golden, load_golden

kms = u.km / u.s

# np.trapz was removed in NumPy 2.0; tumult carries the same shim (tumult/convolutions.py).
trapz = getattr(np, "trapz", None) or np.trapezoid


def test_matches_golden():
    assert_matches_golden(convolution_cases(), load_golden("convolutions"))


def test_kernels_are_normalised():
    from tumult.convolutions import combined_kernel

    px, py = combined_kernel(rd=2e-6, vdop=3.0 * kms, vrot=2.0 * kms, zeta_rt=1.5 * kms)
    area = trapz(np.asarray(py, float), x=np.asarray(px, float))
    assert area == pytest.approx(1.0, rel=1e-6)


def test_convolution_conserves_equivalent_width():
    """Broadening redistributes flux; it must conserve the line's EW (physics invariant)."""
    from tumult.spectrum import Spectrum

    x, y = gauss_line(n=4001, span=4.0)
    sp = Spectrum(x, y)
    before = sp.get_weq(x0=5000.0 * u.AA, profile="trapz", return_profile=True).maa
    conv = sp.crop(vdop=2.0, vrot=1.5, zeta_rt=1.0)
    after = conv.get_weq(x0=5000.0 * u.AA, profile="trapz", return_profile=True).maa
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


def test_custom_kernel_fixed_and_normalised():
    """Phase 1 fixed two custom-kernel bugs: combined_kernel passed 3 args to the
    2-arg custom_kernel (TypeError), and the kernel was never renormalised on the
    convolution grid (leaving a ~1000x composite-unit scale error). The feature was
    unusable before, so there is no golden to preserve - pin the corrected physics:
    a normalised kernel whose width matches the requested velocity profile."""
    from tumult.convolutions import combined_kernel

    vgrid = np.linspace(-10, 10, 201)        # +-5 sigma: truncation bias < 2e-5
    xc = vgrid * kms
    yc = np.exp(-0.5 * (vgrid / 2.0) ** 2)   # sigma = 2 km/s in velocity space
    px, py = combined_kernel(rd=2e-7, custom=(xc, yc))
    px = np.asarray(px, float)
    py = np.asarray(py, float)
    assert trapz(py, x=px) == pytest.approx(1.0, rel=1e-6)
    # second moment of the kernel = (sigma/c)^2 for the Gaussian input
    var = trapz(py * px**2, x=px)
    sigma_over_c = 2.0 / 299792.458
    assert np.sqrt(var) == pytest.approx(sigma_over_c, rel=1e-2)


# --- optimal-filter preparation and resolution bounds (moved out of the filter tool) -------

def test_brault_prep_removes_the_continuum_and_tapers():
    """The filter models a signal about zero, and a finite spectrum must not ring at its ends."""
    import numpy as np
    from tumult import Spectrum
    from tumult.utils import brault_prep

    x = np.linspace(5000.0, 5010.0, 2001)
    spec = Spectrum(x, 1.0 - 0.3 * np.exp(-0.5 * ((x - 5005.0) / 0.05) ** 2))
    xx, y0, w, d = brault_prep(spec)
    assert np.allclose(xx, x)
    assert y0.max() <= 0.0 + 1e-12          # continuum at 1 removed -> absorption is negative
    assert w[0] == 0.0 and w[-1] == 0.0     # tapered to zero at both ends
    assert d > 0                            # ln-lambda spacing


def test_resolution_is_clamped_to_nyquist_and_the_physical_bounds():
    import numpy as np
    from tumult.utils import R_MAX, R_MIN, clamp_resolution, nyquist_resolution

    d = 1e-6
    nyq = nyquist_resolution(d)
    assert nyq == 0.5 / d
    # a request above Nyquist cannot sharpen anything - there is no information up there
    assert clamp_resolution(1e12, d) == pytest.approx(min(nyq, R_MAX))
    assert clamp_resolution(1.0, d) == pytest.approx(R_MIN)
    assert clamp_resolution(5e3, d) == pytest.approx(5e3)


def test_nyquist_resolution_ignores_sign():
    from tumult.utils import nyquist_resolution
    assert nyquist_resolution(-1e-6) == nyquist_resolution(1e-6)
