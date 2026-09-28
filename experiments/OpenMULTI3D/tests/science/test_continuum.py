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


@pytest.mark.parametrize('coeff', [[1.13], [0.04, 1.13], [0.035, -0.02, 1.08]])
def test_joint_polynomial_continuum_and_velocity_fit(coeff):
    """Recover continuum shape and Doppler shift without changing either input."""
    xx = np.linspace(4998, 5002, 1201)
    model = Spectrum(xx, 1 - .4 * np.exp(-.5 * ((xx - 5000) / .1) ** 2))
    target = model.crop(window=(4999, 5001), vshift=.8, continuum=coeff)
    original = target.ydata.copy()
    best, fitted = model.fit(target, fit=['vshift', 'continuum'],
                             p0={'continuum': [0.] * (len(coeff) - 1) + [1.]},
                             xtol={'vshift': 1e-4, 'continuum': 1e-6})
    assert best['continuum'] == pytest.approx(coeff, abs=2e-5)
    assert best['vshift'] == pytest.approx(.8, abs=.002)
    assert np.max(np.abs(fitted.ydata.value - target.ydata.value)) < 1e-5
    np.testing.assert_array_equal(target.ydata, original)


def test_default_continuum_and_explicit_coefficient_bounds():
    xx = np.linspace(5000, 5002, 101)
    model = Spectrum(xx, np.ones_like(xx))
    target = model.crop(continuum=[.03, 1.07])
    best, _ = model.fit(target, fit=['continuum'])
    assert best['continuum'] == pytest.approx([.03, 1.07], abs=3e-5)
    constrained, _ = model.fit(target, fit=['continuum'],
                               bounds={'continuum': ([-.01, .95], [.01, 1.05])})
    assert constrained['continuum'] == pytest.approx([.01, 1.05], abs=3e-5)


def test_polynomial_domain_is_reusable_and_nodes_follow_pchip():
    xx = np.linspace(4998, 5002, 101)
    spec = Spectrum(xx, np.ones_like(xx))
    domain = (4999, 5001)
    full = spec.crop(continuum=[.1, 1], continuum_window=domain)
    narrow = spec.crop(window=domain, continuum=[.1, 1])
    np.testing.assert_allclose(full.crop(window=domain).ydata, narrow.ydata)
    nodes = [(4999, .9), (5000, 1.1), (5001, 1)]
    out = spec.crop(continuum=nodes)
    np.testing.assert_allclose(out.ydata.value, cont.evaluate(*np.asarray(nodes).T, xx))


def _local_measurement():
    from tumult.abundance import AbundanceResult, LineResult
    from tumult.linemask import LineMask
    xx = np.linspace(4999, 5001, 1001)
    model = Spectrum(xx, 1 - .3 * np.exp(-.5 * ((xx - 5000) / .1) ** 2))
    target = model.crop(continuum=[.015, 1.03])
    ln = LineResult(5000, 'Fe I', 1., -1.,
                    mask=LineMask(lmin=4999.5, lmax=5000.5, smin=4999, smax=5001))
    res = AbundanceResult(symbol='Fe', summary={}, lines=[ln], profile='trapz', snr=100, cont_pct=.5)
    return model, target, ln, res


def test_local_continuum_changes_ew_and_model_but_never_target_or_mask():
    from copy import deepcopy
    from tumult import abundance as ab
    model, target, ln, res = _local_measurement()
    original, mask = target.ydata.copy(), ln.mask.to_dict()
    independent = deepcopy(res)
    ab.set_local_continuum(res, ln, target, enabled=True)
    assert np.asarray(ln.continuum)[:, 0] == pytest.approx(ln.window)
    assert ln.ew == pytest.approx(ab.measure_ew(model, ln.fit_window).maa, rel=.002)
    first_ew = ln.ew
    nodes = ((4999., 1.025), (5001., 1.055))
    ab.set_local_continuum(res, ln, target, nodes=nodes)
    assert ln.ew > first_ew
    ln.fit_spec = model
    np.testing.assert_allclose(ab.displayed_model(ln).ydata.value,
                               model.ydata.value * cont.model_continuum(nodes, model.values[0]))
    assert independent.lines[0].renorm is False
    assert ln.mask.to_dict() == mask
    np.testing.assert_array_equal(target.ydata, original)
    restored = ab.AbundanceResult.from_dict(res.to_dict())
    assert restored.lines[0].continuum == nodes
    assert restored.measurement().lines[0].continuum == nodes
    ab.set_local_continuum(res, ln, target, enabled=False)
    assert ln.ew == pytest.approx(ab.measure_ew(target, ln.fit_window).maa)
    np.testing.assert_array_equal(target.ydata, original)


def test_edited_ew_requests_only_missing_nearest_rungs():
    from tumult import abundance as ab
    model, target, ln, res = _local_measurement()
    xx = model.values[0]
    def spectrum(a):
        return Spectrum(xx, 1 - .2 * 10**(.4*a) * np.exp(-.5*((xx-5000)/.1)**2))
    ln.method = 'cog'
    ln.cog_specs = {a: spectrum(a) for a in [-1., 0., 1.]}
    ln.cog_step = .25
    ln.ew = ab.measure_ew(spectrum(.62), ln.fit_window).maa
    ln.ew_err = 1.
    original = dict(ln.cog_specs)
    ab.refresh_line_cog(ln, 'trapz')
    assert ln.missing_rungs == (.5, .75)
    assert 'missing_rungs' in ln.flags and np.isnan(ln.a) and ln.include
    calls = []
    def synthesize(values):
        calls.extend(values)
        return [[spectrum(a)] for a in values]
    assert ab.complete_cog_rungs([ln], synthesize, profile='trapz')
    assert calls == [.5, .75]
    assert ln.a == pytest.approx(.62, abs=1e-6)
    assert ln.missing_rungs == () and 'missing_rungs' not in ln.flags
    assert all(ln.cog_specs[a] is sp for a, sp in original.items())
    # Moving back into an already cached cell needs no synthesis.
    ln.ew = ab.measure_ew(spectrum(.55), ln.fit_window).maa
    ab.refresh_line_cog(ln, 'trapz')
    assert ln.missing_rungs == () and ln.a == pytest.approx(.55, abs=1e-6)


def test_local_edit_without_rungs_invalidates_a_direct_fit():
    from tumult import abundance as ab
    model, target, ln, res = _local_measurement()
    ln.a, ln.drewda, ln.sigma_a, ln.fit_spec = .1, .4, .02, model
    ab.set_local_continuum(res, ln, target, enabled=True)
    assert 'refit_required' in ln.flags
    assert np.isnan(ln.a) and np.isnan(ln.sigma_a) and np.isnan(ln.drewda)
    assert ln.fit_spec is not None
    np.testing.assert_allclose(ab.displayed_model(ln).ydata.value,
                               model.ydata.value * cont.model_continuum(
                                   ln.continuum, model.values[0]))
