"""Pin the continuum placement extracted out of the canvas and the normalise tool.

The properties that matter are the ones a wrong continuum would break downstream: the level
must be robust to spikes, the spline must not overshoot between knots (an overshoot puts flux
above the continuum and turns a line into emission), it must hold flat rather than diverge
past the outermost knot, and the knot list handed to pchip must be strictly increasing even
when local windows drop colinear knots on top of the global ones.
"""
import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d

from tumult import continuum as cont
from tumult import Spectrum
from tumult.utils import pchip_mono


def test_flat_level_ignores_spikes_and_nans():
    y = np.r_[np.ones(200), 50.0, np.nan]
    assert cont.flat_level(y) == pytest.approx(1.0, abs=1e-9)


def test_flat_continuum_broadcasts():
    c = cont.flat_continuum(np.ones(7) * 3.0)
    assert c.shape == (7,)
    assert np.allclose(c, 3.0)


def test_merge_knots_defaults_to_a_flat_span():
    kx, ky = cont.merge_knots(np.linspace(0, 10, 50), np.ones(50))
    assert (kx[0], kx[-1]) == pytest.approx((0.0, 10.0))
    assert np.allclose(ky, 1.0)


def test_merge_knots_lets_a_window_govern_its_interior():
    """A base knot strictly inside a renorm window is dropped for the window's own knots."""
    kx, ky = cont.merge_knots(
        [0, 10], [1, 1], base_knots=([0.0, 5.0, 10.0], [1.0, 1.0, 1.0]),
        renorm_windows=[(4.0, 6.0, [(4.5, 2.0), (5.5, 2.0)])])
    assert 5.0 not in list(kx)
    assert list(kx) == pytest.approx([0.0, 4.5, 5.5, 10.0])
    assert list(ky) == pytest.approx([1.0, 2.0, 2.0, 1.0])


def test_merge_knots_is_strictly_increasing_on_collisions():
    """pchip_mono needs increasing x; colinear window knots readily land on a base knot."""
    kx, _ky = cont.merge_knots(
        [0, 10], [1, 1], base_knots=([0.0, 10.0], [1.0, 1.0]),
        renorm_windows=[(1.0, 9.0, [(5.0, 2.0), (5.0, 2.0), (5.0, 2.0)])])
    assert np.all(np.diff(kx) > 0)


def test_evaluate_holds_flat_beyond_the_knots():
    """Past the outermost knot the continuum must not run away on the cubic tail."""
    kx, ky = np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 1.5])
    out = cont.evaluate(kx, ky, np.array([-100.0, 1.0, 3.0, 100.0]))
    assert out[0] == pytest.approx(out[1])
    assert out[-1] == pytest.approx(out[-2])


def test_evaluate_does_not_overshoot_between_knots():
    kx = np.array([0.0, 1.0, 2.0, 3.0])
    ky = np.array([1.0, 1.0, 2.0, 2.0])
    y = cont.evaluate(kx, ky, np.linspace(0, 3, 400))
    assert y.min() >= ky.min() - 1e-9
    assert y.max() <= ky.max() + 1e-9


def test_evaluate_matches_a_clamped_pchip():
    """Equivalence with the expression this replaced in the canvas."""
    kx, ky = np.array([1.0, 2.5, 4.0]), np.array([1.0, 1.4, 0.9])
    x = np.linspace(-1, 6, 200)
    assert np.allclose(cont.evaluate(kx, ky, x),
                       pchip_mono(kx, ky)(np.clip(x, kx.min(), kx.max())))


def test_sanitise_blanks_zero_and_nonfinite():
    out = cont.sanitise([1.0, 0.0, np.inf, -np.inf, np.nan, 2.0])
    assert np.isnan(out[1:5]).all()
    assert list(out[[0, 5]]) == [1.0, 2.0]


def _spiky():
    x = np.linspace(5000.0, 5010.0, 2000)
    y = 1.0 - 0.4 * np.exp(-0.5 * ((x - 5005.0) / 0.05) ** 2)
    y[::97] += 0.3                       # spikes above the continuum
    return Spectrum(x, y)


def test_reject_outliers_drops_points_above_the_running_mean():
    spec = _spiky()
    out = cont.reject_outliers(spec, window=0.5, percentile=80.0)
    assert len(out.values[0]) < len(spec.values[0])


def test_reject_outliers_matches_the_inline_form():
    """Equivalence with the expression this replaced in the normalise tool."""
    spec, win, pct = _spiky(), 0.5, 80.0
    x, y = spec.values
    dx = abs(float(np.mean(np.diff(x))))
    nstep = max(1, int(round(win / dx)))
    mask = (y - uniform_filter1d(y, nstep, mode="nearest")) <= np.percentile(
        y - uniform_filter1d(y, nstep, mode="nearest"), pct)
    expected = x[mask]
    assert np.allclose(cont.reject_outliers(spec, window=win, percentile=pct).values[0],
                       expected)


def test_reject_outliers_is_a_noop_without_settings():
    spec = _spiky()
    assert len(cont.reject_outliers(spec).values[0]) == len(spec.values[0])
    assert len(cont.reject_outliers(spec, window=0, percentile=0).values[0]) == len(spec.values[0])


def test_reject_outliers_keeps_the_spectrum_when_nothing_survives():
    """A percentile that would leave under two points must not empty the spectrum."""
    spec = _spiky()
    out = cont.reject_outliers(spec, window=0.5, percentile=1e-9)
    assert len(out.values[0]) == len(spec.values[0])
