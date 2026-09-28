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


def test_local_continuum_edit_refits_a_profile_fit_on_its_rungs():
    """Editing a profile-fitted line's local continuum fits it again on the rungs the fit
    kept, without synthesis, and switching the continuum off restores the first solution."""
    nodes = ((4999.2, 1.02), (5000.8, 1.06))
    line = LineResult(lam0=5000.0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=90.0, ew_err=3.0,
                      fwhm=0.19, mask=_mask())
    target = _spectrum(0.35 * 10 ** (0.5 * 0.23)).crop(continuum=list(nodes))
    result = ab.AbundanceResult(symbol="Fe", summary={}, lines=[line], profile="gauss",
                                snr=100, cont_pct=0.0)
    calls = []

    def synthesize(offsets):
        calls.extend(offsets)
        return [[_spectrum(0.35 * 10 ** (0.5 * a))] for a in offsets]

    ab.fit_profiles([line], synthesize, target, step=0.25, limit=1.0, max_rounds=8,
                    profile="gauss")
    first, computed = line.a, len(calls)
    assert set(line.cog_specs) == set(calls)

    ab.set_local_continuum(result, line, target, nodes=nodes, enabled=True)
    assert len(calls) == computed
    assert not {"refit_required", "missing_rungs", "out_of_range"} & set(line.flags)
    assert line.a == pytest.approx(0.23, abs=0.02)
    assert np.isfinite(line.drewda) and np.isfinite(line.sigma_a)
    assert result.summary["a_mean"] == pytest.approx(line.a)

    ab.set_local_continuum(result, line, target, enabled=False)
    assert len(calls) == computed
    assert line.a == pytest.approx(first, abs=0.01)


@pytest.mark.parametrize("method", ["fit", "cog"])
def test_local_continuum_divides_the_synthesis_by_its_own_pseudo_continuum(method):
    """With a local continuum on, each rung is divided by its own straight continuum, so a
    tilted synthetic pseudo-continuum does not bias the abundance; with it off, it does."""
    nodes = ((4999.2, 1.02), (5000.8, 1.06))
    x = np.arange(4999.0, 5001.0, 0.01)
    tilt = 0.94 + 0.03 * (x - 4999.0)       # 3% below unity at the line

    def synthesize(offsets):
        return [[_spectrum(0.35 * 10 ** (0.5 * a)) * tilt] for a in offsets]

    target = _spectrum(0.35 * 10 ** (0.5 * 0.23)).crop(continuum=list(nodes))
    line = LineResult(lam0=5000.0, ion="Fe I", e_low=1.0, loggf=-1.0, fwhm=0.19,
                      mask=_mask(), renorm=True, continuum=nodes)
    result = ab.AbundanceResult(symbol="Fe", summary={}, lines=[line], profile="gauss",
                                snr=100, cont_pct=0.0)
    ab.set_local_continuum(result, line, target)
    if method == "fit":
        ab.fit_profiles([line], synthesize, target, step=0.25, limit=1.0, max_rounds=8,
                        profile="gauss")
        tolerance = 0.02
    else:
        ab.adaptive_cog([line], synthesize, step=0.25, limit=1.0, max_rounds=8,
                        profile="gauss")
        tolerance = 0.03
    assert line.a == pytest.approx(0.23, abs=tolerance)
    assert line.fit_spec.ydata.value.max() == pytest.approx(1.0, abs=0.01)

    ab.set_local_continuum(result, line, target, enabled=False)
    assert abs(line.a - 0.23) > tolerance or np.isnan(line.a)
