"""pchip_mono with Quantity inputs.

The shape-preserving end-segment requirement moved to
`tests/science/test_interpolation.py` as part of requirement S11 of
`reports/test-requirements.md`. The unit handling below has no approved requirement of
its own yet.
"""
import numpy as np
import astropy.units as u

from tumult.utils import pchip_mono


def test_quantity_inputs_still_interpolate():
    x = np.linspace(5000.0, 5001.0, 6) * u.AA
    y = np.exp(np.linspace(0.0, 1.0, 6)) * u.dimensionless_unscaled
    out = pchip_mono(x, y)(x[:2])
    assert np.allclose(np.asarray(out), np.asarray(y[:2]))
