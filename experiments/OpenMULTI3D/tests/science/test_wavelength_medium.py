"""Air and vacuum wavelengths: the round trip, and the size of the refraction offset.

Requirement S6 of `reports/test-requirements.md`.

Assumptions:

* wavelengths are in Angstrom throughout, and the conversion is the one for *standard
  air*: dry air at 15 degrees Celsius and 101.325 kPa with the usual 0.045% CO2. The
  dispersion formula carries those conditions in its coefficients, so a spectrum taken
  at another temperature or pressure is not converted exactly by it. Every user-facing
  wavelength in OpenMULTI3D is on that air scale and every native grid is in vacuum;
  see *wavelength medium* in `docs/GLOSSARY.md`;
* the span covered is 3000 to 25000 Angstrom, the optical and near infrared the code is
  used in. The formula is a fit and is not claimed outside it.
"""
import numpy as np
import pytest

# The Ba II resonance line, whose vacuum wavelength is a convenient anchor for the size
# of the offset in the blue-green.
BA_II_VACUUM_AA = 4554.033

GRID_AA = np.linspace(3000.0, 25000.0, 2001)


def test_air_to_vacuum_round_trips_over_the_optical_and_near_infrared():
    """Converting a wavelength there and back returns it.

    The vacuum-to-air direction divides by the refractive index evaluated at the vacuum
    wavelength; the inverse has no closed form and is reached by a fixed number of
    fixed-point steps. The round trip therefore tests the convergence of that iteration.

    Tolerance: atol 1e-8 Angstrom, absolute because the quantity being compared is a
    difference that should be zero, and a relative bound near 25000 Angstrom would be a
    far weaker statement than near 3000. The bound is set by the iteration and not by
    floating point: n - 1 is about 3e-4, so each fixed-point step reduces the residual
    by that factor and three steps from a starting error of order 1 Angstrom leave about
    1e-11 Angstrom. The bound leaves three decades of headroom over that, and would
    still fail if a step were dropped.
    """
    from tumult.physics import air2vac, vac2air

    np.testing.assert_allclose(air2vac(vac2air(GRID_AA)), GRID_AA, atol=1e-8, rtol=0)
    np.testing.assert_allclose(vac2air(air2vac(GRID_AA)), GRID_AA, atol=1e-8, rtol=0)


def test_vacuum_wavelengths_exceed_air_wavelengths_everywhere():
    """Air refracts, so n > 1 and lambda_air = lambda_vac / n is the shorter of the two.

    This is the sign of the effect, and it is the defect most likely to go unnoticed:
    an inverted conversion is a shift of one to four Angstrom, which looks like a
    velocity rather than like a bug. The offset also grows with wavelength across the
    span, since n approaches 1 towards the infrared but lambda (n - 1) does not.
    """
    from tumult.physics import vac2air

    air = vac2air(GRID_AA)
    assert np.all(air < GRID_AA)
    offsets = GRID_AA - air
    assert np.all(np.diff(offsets) > 0)


def test_the_offset_near_4554_angstrom_matches_the_literature_value():
    """lambda_vac - lambda_air = 1.28 Angstrom at 4554 Angstrom.

    Derivation, independent of the code: with sigma = 1/lambda_vac in inverse microns,
    the standard-air dispersion of Birch and Downs (1994) gives
    (n - 1) * 1e8 = 8342.54 + 2406147/(130 - sigma^2) + 15998/(38.9 - sigma^2).
    At lambda_vac = 4554.033 AA, sigma = 2.19587 and sigma^2 = 4.82184, so
      2406147/125.178 = 19221.2 and 15998/34.078 = 469.5,
    giving (n - 1) = 2.80332e-4. The offset is lambda_vac (1 - 1/n) = lambda_air (n - 1)
    = 4552.76 * 2.80332e-4 = 1.2763 Angstrom.

    Tolerance: atol 0.05 Angstrom, absolute because the comparison is against a value
    quoted to two decimals. The band is about 4% of the offset: wide enough to accept
    any of the standard-air formulations in use (Edlen 1966, Birch and Downs 1994,
    Ciddor 1996 differ by a few times 1e-8 in n, that is under 1e-3 Angstrom here) and
    narrow enough to reject a formula evaluated at the wrong argument or with a missing
    dispersion term, which shift the offset by several percent.
    """
    from tumult.physics import vac2air

    air = float(vac2air(BA_II_VACUUM_AA))
    assert BA_II_VACUUM_AA - air == pytest.approx(1.28, abs=0.05)
