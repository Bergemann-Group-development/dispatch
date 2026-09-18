"""pchip_mono must be shape preserving at the ends as well as in the interior.

A three-rung NLTE curve of growth with a saturated first rung (REW spacing 0.11 then
0.03 dex per dex of abundance) made the unguarded three-point endpoint slope change sign,
and the inverted abundance dipped below the first rung before rising through the data.
"""
import numpy as np
import astropy.units as u

from tumult.utils import pchip_mono


REW = np.array([-4.812, -4.700, -4.672])
DA = np.array([-1.0, 0.0, 1.0])


def test_end_segment_stays_between_its_rungs():
    xs = np.linspace(REW[0], REW[1], 50)
    a = pchip_mono(REW, DA)(xs)
    assert a.min() >= DA[0] - 1e-12
    assert a.max() <= DA[1] + 1e-12
    assert np.all(np.diff(a) >= 0)


def test_quantity_inputs_still_interpolate():
    x = np.linspace(5000.0, 5001.0, 6) * u.AA
    y = np.exp(np.linspace(0.0, 1.0, 6)) * u.dimensionless_unscaled
    out = pchip_mono(x, y)(x[:2])
    assert np.allclose(np.asarray(out), np.asarray(y[:2]))
