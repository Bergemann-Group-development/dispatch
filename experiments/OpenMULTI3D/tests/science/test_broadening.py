"""Broadening and shifting: kernel normalisation, width conservation, the Doppler shift.

Requirements S7, S8, S9 and S10 of `reports/test-requirements.md`.

Shared assumptions:

* a kernel's abscissa is the dimensionless `delta lambda / lambda`, which for a
  non-relativistic velocity is `v/c`. The kernels are therefore built and normalised in
  velocity space and applied on a logarithmic wavelength grid, where a velocity is a
  constant shift;
* `rd` is the step of that logarithmic grid, also dimensionless;
* spectra are continuum normalised with a continuum of exactly 1.0, and wavelengths are
  in Angstrom on the air scale (the user-facing medium; see *wavelength medium* in
  `docs/GLOSSARY.md`). Broadening and shifting act on the axis as given, so no
  conversion takes place and the results do not depend on the medium;
* the velocity broadening applied is small compared with the span of the spectrum, which
  is what keeps the flux a kernel moves inside the window rather than off its edge.
"""
import astropy.units as u
import numpy as np
import pytest

from tests.support.spectra import (LINE_CENTRE, LINE_LAM0, flux_centroid, gaussian_line,
                                   gaussian_line_equivalent_width_maa)

KMS = u.km / u.s
C_KMS = 299792.458          # speed of light, the value astropy's constant carries
RD = 2e-6                   # the ln-lambda step the kernels below are built on

trapezoid = getattr(np, "trapezoid", None) or np.trapz


def test_composed_kernel_integrates_to_one_on_its_own_abscissa():
    """A broadening kernel redistributes flux, so its area is 1 by construction.

    The composition is rotation (2 km/s), a Gaussian Doppler term (3 km/s) and
    radial-tangential macroturbulence (1.5 km/s). Each is normalised separately and the
    three are then convolved together, which is where the invariant can be lost: the
    convolution of two unit-area kernels has unit area only if the discrete sum carries
    the grid step, and a missing step factor scales the result by 1/rd.

    Tolerance: rel 1e-6, the trapezoidal error of integrating the composed kernel on the
    same grid it is sampled on. That grid resolves the narrowest component with many
    points, so the quadrature error is far below the bound; a lost step factor is a
    factor of 500000 away from it.
    """
    from tumult.convolutions import combined_kernel

    px, py = combined_kernel(rd=RD, vdop=3.0 * KMS, vrot=2.0 * KMS, zeta_rt=1.5 * KMS)
    area = trapezoid(np.asarray(py, float), x=np.asarray(px, float))
    assert area == pytest.approx(1.0, rel=1e-6)


def test_custom_kernel_integrates_to_one_on_the_convolution_grid():
    """A caller-supplied kernel is resampled onto the convolution grid and renormalised.

    This is a different construction from the analytic kernels above: the shape arrives
    as (velocity, value) samples of the caller's own choosing, is interpolated onto the
    grid and must come back with unit area there. Interpolation alone does not preserve
    an area, and the input carries a velocity unit that has to leave the kernel, so the
    renormalisation is load-bearing rather than redundant.

    The input is a Gaussian of 2 km/s standard deviation sampled over +-5 sigma, where
    truncation removes 6e-7 of its area.
    """
    from tumult.convolutions import combined_kernel

    velocity = np.linspace(-10.0, 10.0, 201)
    shape = np.exp(-0.5 * (velocity / 2.0) ** 2)
    px, py = combined_kernel(rd=2e-7, custom=(velocity * KMS, shape))
    area = trapezoid(np.asarray(py, float), x=np.asarray(px, float))
    assert area == pytest.approx(1.0, rel=1e-6)


def test_gaussian_velocity_kernel_has_the_requested_second_moment():
    """The second moment of the Doppler kernel is (vdop/c)^2 / 2.

    Derivation: the Gaussian Doppler kernel is written with the 1/e half width
    b = vdop/c, that is exp(-(x/b)^2), not with a standard deviation. Its normalised
    second moment is
      integral x^2 exp(-(x/b)^2) dx / integral exp(-(x/b)^2) dx = b^2/2,
    so the root mean square abscissa is b/sqrt(2). For vdop = 3 km/s this is
    3/299792.458/sqrt(2) = 7.07596e-6 in delta lambda / lambda, that is 2.12132 km/s
    expressed as a velocity.

    The convention matters more than the number: reading vdop as a standard deviation
    instead of a 1/e width is a factor sqrt(2) in every line width the code produces.

    Tolerance: rel 1e-2, set by the truncation of the kernel at its finite support and
    by the trapezoidal error of the moment on a grid of step rd. The bound is far wider
    than either, and is chosen so that it still separates the two conventions, which
    differ by 41%.
    """
    from tumult.convolutions import combined_kernel

    px, py = combined_kernel(rd=RD, vdop=3.0 * KMS)
    px = np.asarray(px, float)
    py = np.asarray(py, float)
    rms = np.sqrt(trapezoid(py * px**2, x=px))
    assert rms == pytest.approx(3.0 / C_KMS / np.sqrt(2.0), rel=1e-2)


def test_convolution_conserves_the_equivalent_width():
    """Broadening moves flux within the profile; the area between line and continuum is kept.

    This is the invariant that makes an equivalent width a measure of how much light a
    line removes rather than of how the star smears it: convolving a normalised spectrum
    with a unit-area kernel leaves the integral of (1 - y) unchanged, provided none of
    the redistributed flux leaves the window.

    Padding: the line sits at 5000.05 Angstrom with sigma 0.11 Angstrom in a window
    reaching 4 Angstrom each side of 5000. The composed kernel (2 km/s Doppler, 1.5 km/s
    rotation, 1 km/s macroturbulence) has a half support of 0.225 Angstrom here, so its
    full width at half maximum is under 0.45 Angstrom. The distance from the far wing of
    the line (four of its own sigma, 0.44 Angstrom) to either edge is about 3.5
    Angstrom, more than five times that full width at half maximum.

    Tolerance: rel 1e-3, covering the flux carried across the window edge and the
    interpolation onto the logarithmic grid and back. The measured change is far under
    it; the bound is set where it would still catch an unnormalised kernel or a
    convolution that truncates the wings.
    """
    from tumult.spectrum import Spectrum

    x, y = gaussian_line(n_points=4001, half_span=4.0)
    spectrum = Spectrum(x, y)
    broadened = spectrum.crop(vdop=2.0, vrot=1.5, zeta_rt=1.0)

    width = broadened.get_weq(x0=LINE_LAM0 * u.AA, profile="trapz",
                              return_profile=True).maa
    assert width == pytest.approx(gaussian_line_equivalent_width_maa(), rel=1e-3)


def test_a_velocity_shift_moves_the_feature_and_leaves_the_grid_alone():
    """lambda' = lambda (1 + v/c), with a positive velocity moving the line redward.

    Derivation: the non-relativistic Doppler relation puts a feature of rest wavelength
    5000.05 Angstrom at 5000.05 * (1 + 10/299792.458) = 5000.21678 Angstrom for a
    receding source at 10 km/s, a displacement of 0.16678 Angstrom.

    The recovered centre is the absorption-weighted centroid of the shifted spectrum,
    which is exact for a symmetric profile on a uniform grid and is computed here
    without the interpolation the shift itself uses.

    Tolerance: half a sample spacing, 0.0025 Angstrom, which is the resolution at which
    a position on this grid is meaningful at all. That is a sixtieth of the displacement
    being tested, so a sign error or a missing (1 + v/c) factor cannot pass.

    The wavelength grid is asserted unchanged with exact equality: the shift resamples
    the flux onto the axis it was given, and a caller relies on being able to compare
    the two spectra point by point.
    """
    from tumult.spectrum import Spectrum

    x, y = gaussian_line()
    spectrum = Spectrum(x, y)
    shifted = spectrum.crop(vshift=10.0)

    np.testing.assert_array_equal(np.asarray(shifted.xdata.value, float), x)

    centre = flux_centroid(np.asarray(shifted.xdata.value, float),
                           np.asarray(shifted.ydata, float))
    spacing = x[1] - x[0]
    assert centre == pytest.approx(LINE_CENTRE * (1.0 + 10.0 / C_KMS),
                                   abs=0.5 * spacing)
    assert centre > flux_centroid(x, y)


# --- resolving power ------------------------------------------------------------------
#
# These are discrete decisions about a requested number, so they are asserted exactly.
# R_MIN and R_MAX are the physical bounds (below R_MIN a filter removes the spectrum
# rather than its noise; above R_MAX there is nothing there but noise) and the Nyquist
# limit of a logarithmic grid of step d is 0.5/d, two samples per resolution element.


def test_nyquist_resolution_is_half_the_inverse_step_and_ignores_the_sign():
    """nyquist = 0.5/|d|. A descending grid samples exactly as finely as an ascending one.

    Exact equality: this is a definition, not a measurement.
    """
    from tumult.utils import nyquist_resolution

    assert nyquist_resolution(1e-6) == 0.5 / 1e-6
    assert nyquist_resolution(-1e-6) == nyquist_resolution(1e-6)


def test_a_request_above_the_nyquist_limit_is_clamped_to_it():
    """Sampling of step 1e-6 carries nothing above R = 500000, so that is what is returned.

    There is no information above the Nyquist limit to keep, so the request is clamped
    rather than honoured. Exact equality, for the same reason as above.
    """
    from tumult.utils import R_MAX, clamp_resolution

    nyquist = 0.5 / 1e-6
    assert nyquist < R_MAX
    assert clamp_resolution(1e12, 1e-6) == nyquist


def test_a_request_above_the_physical_bound_is_clamped_to_it():
    """Where the sampling is finer than the physics, R_MAX is the binding limit.

    With a step of 1e-8 the Nyquist limit is 5e7, well above R_MAX, so a very large
    request comes back as R_MAX. This separates the two clamps, which the case above
    cannot: with only one grid step both would be satisfied by whichever limit happened
    to be applied.
    """
    from tumult.utils import R_MAX, clamp_resolution

    assert 0.5 / 1e-8 > R_MAX
    assert clamp_resolution(1e12, 1e-8) == R_MAX


def test_a_request_below_the_physical_bound_is_clamped_up_and_one_inside_is_kept():
    """R_MIN is the floor, and a request between the bounds passes through unchanged."""
    from tumult.utils import R_MIN, clamp_resolution

    assert clamp_resolution(1.0, 1e-6) == R_MIN
    assert clamp_resolution(5e3, 1e-6) == 5e3


def test_the_clamp_does_not_depend_on_the_sign_of_the_step():
    """A descending wavelength axis must not silently lose the Nyquist limit.

    The step of a descending grid is negative, and the limit is a property of the
    spacing rather than of the direction. Exact equality between the two.
    """
    from tumult.utils import clamp_resolution

    assert clamp_resolution(1e12, -1e-6) == clamp_resolution(1e12, 1e-6)
    assert clamp_resolution(5e3, -1e-6) == clamp_resolution(5e3, 1e-6)
