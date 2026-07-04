"""Pin the air/vacuum wavelength conversion (single-source Ciddor-style formula in
auxfuncs; moves to tumult/physics.py in Phase 1)."""
import numpy as np
import pytest

from cases import AIRVAC_LINES, airvac_cases
from conftest import assert_matches_golden, load_golden


def test_matches_golden():
    assert_matches_golden(airvac_cases(), load_golden("airvac"))


def test_round_trip():
    """air2vac(vac2air(x)) must return x to < 1e-8 AA across the optical/NIR."""
    from tumult.auxfuncs import air2vac, vac2air

    grid = np.linspace(3000.0, 25000.0, 2001)
    np.testing.assert_allclose(air2vac(vac2air(grid)), grid, atol=1e-8, rtol=0)
    np.testing.assert_allclose(vac2air(air2vac(grid)), grid, atol=1e-8, rtol=0)


def test_refraction_sign_and_magnitude():
    """n_air > 1: vacuum wavelengths are longer; the shift is ~1.3 AA near 4554 AA."""
    from tumult.auxfuncs import vac2air

    air = vac2air(AIRVAC_LINES)
    assert np.all(air < AIRVAC_LINES)
    shift = AIRVAC_LINES - air
    assert 1.0 < shift[1] < 1.6  # Ba II 4554: literature air-vac offset ~1.28 AA
