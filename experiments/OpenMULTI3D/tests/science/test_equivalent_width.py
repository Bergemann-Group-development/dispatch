"""Equivalent width: the analytic width of a Gaussian line, its units, its error budget.

Requirements S1, S2 and S3 of `reports/test-requirements.md`.

Shared assumptions, which hold for every case in this module:

* the flux is already continuum normalised, so the continuum is exactly 1.0 and no
  continuum placement enters the measurement;
* wavelengths are in Angstrom on the air scale (the user-facing medium; see
  *wavelength medium* in `docs/GLOSSARY.md`). The width and the reduced width are
  read off the same axis the line was built on, so neither depends on the medium;
* widths are in milli-Angstrom wherever the attribute name says `maa`, and the
  reduced equivalent width `log10(W / lambda)` is dimensionless;
* the line is a single unblended feature in the weak-line regime: no saturation, no
  second component, and a window of 18 sigma each side at 22 samples per sigma.
"""
import astropy.units as u
import numpy as np
import pytest

from tests.support.spectra import LINE_LAM0, gaussian_line, gaussian_line_equivalent_width_maa

# Every profile the package registers. A profile absent from this list would never be
# measured here, so the set is asserted rather than assumed: S1 is a statement about
# *every* registered profile.
FITTED_PROFILES = ("gauss", "voigt", "asym_gauss")
INTEGRATED_PROFILES = ("trapz", "pchip")


def test_registry_is_exactly_the_profiles_measured_here():
    """The parameterisation below covers the whole registry, so S1 means what it says."""
    from tumult.profile import profile_names

    assert set(profile_names()) == set(FITTED_PROFILES) | set(INTEGRATED_PROFILES)


@pytest.mark.parametrize("profile", FITTED_PROFILES)
def test_fitted_profile_recovers_the_analytic_gaussian_width(profile):
    """W = A sqrt(2 pi) sigma for a Gaussian of depth A and width sigma.

    Derivation: the line is 1 - A exp(-(x - x0)^2 / 2 sigma^2) on a unit continuum, so
    W = integral (1 - y) dx = A sigma sqrt(2 pi). With A = 0.42 and sigma = 0.11 AA
    that is 0.42 * 2.5066283 * 0.11 = 0.11580623 AA = 115.80623 mAA.

    Tolerance: rel 1e-6. The data are noiseless and the three fitted shapes all contain
    the true Gaussian, so the least-squares minimum sits at the exact parameters and the
    only error left is the fitter's own convergence threshold. A bound this tight would
    catch a wrong normalisation constant (sqrt(2 pi) against 2, a 25% error) many orders
    over, and also a Voigt or asymmetric-Gaussian width formula that disagrees with the
    Gaussian limit its extra parameter collapses to.
    """
    from tumult.utils import get_weq

    x, y = gaussian_line()
    measured = get_weq(x, y, x0=LINE_LAM0, profile=profile, return_profile=True)
    assert measured.maa == pytest.approx(gaussian_line_equivalent_width_maa(), rel=1e-6)


@pytest.mark.parametrize("profile", INTEGRATED_PROFILES)
def test_integrated_profile_recovers_the_analytic_gaussian_width(profile):
    """Same analytic width, reached by integrating the sampled line instead of fitting it.

    Tolerance: rel 1e-3, the sum of two errors of the integration. The trapezoidal rule
    at 22 samples per sigma is far better than its O(h^2) worst case on this integrand,
    because a Gaussian and its derivatives vanish at both ends of the window; the
    truncated wings beyond 18 sigma contribute exp(-162) of the depth. Both are far
    under the bound, which is set instead at the level that still separates a correct
    width from a wrong quadrature weight or a missing first or last sample.
    """
    from tumult.utils import get_weq

    x, y = gaussian_line()
    measured = get_weq(x, y, x0=LINE_LAM0, profile=profile, return_profile=True)
    assert measured.maa == pytest.approx(gaussian_line_equivalent_width_maa(), rel=1e-3)


def test_width_converts_consistently_between_wavelength_units():
    """A width is one physical quantity: AA, nm and mAA views differ by the unit factor.

    Derivation: 1 AA = 0.1 nm = 1000 mAA, so W[mAA] = 1000 W[AA] = 10000 W[nm].

    Tolerance: rel 1e-12. Each conversion is a single scalar multiply on a value the
    object already holds, so anything above rounding is a wrong factor.
    """
    from tumult.utils import get_weq

    x, y = gaussian_line()
    measured = get_weq(x, y, x0=LINE_LAM0, profile="trapz", return_profile=True)

    in_aa = measured.width(u.AA)
    assert measured.width(u.mAA) == pytest.approx(in_aa * 1e3, rel=1e-12)
    assert measured.width(u.nm) == pytest.approx(in_aa * 0.1, rel=1e-12)
    assert measured.maa == pytest.approx(measured.width(u.mAA), rel=1e-12)


def test_reduced_equivalent_width_is_log10_of_the_width_over_the_wavelength():
    """REW = log10(W / lambda) with W in mAA and lambda in AA.

    Derivation: 50 mAA at 5000 AA is W/lambda = 0.050/5000 = 1e-5, so REW = -5 exactly.
    A second anchor: 500 mAA at 5000 AA gives 1e-4 and REW = -4.

    Tolerance: rel 1e-12 on a value that is exact in binary except for log10's own
    rounding. The pairing of units is the point: mAA over AA needs the factor 1e-3, and
    omitting it would put the answer three decades out.
    """
    from tumult.profile import reduced_ew

    assert reduced_ew(50.0, 5000.0) == pytest.approx(-5.0, rel=1e-12)
    assert reduced_ew(500.0, 5000.0) == pytest.approx(-4.0, rel=1e-12)


def test_reduced_equivalent_width_guards_non_positive_inputs():
    """A logarithm has no value at or below zero, so the guard returns -inf exactly.

    Exact equality is right here: these are discrete decisions of the guard, not
    measurements. A line measured as zero or negative width (a filled-in or emission
    feature) must read as infinitely weak rather than raise or return NaN, since the
    caller ranks lines by this number.
    """
    from tumult.profile import reduced_ew

    assert reduced_ew(0.0, 5000.0) == -np.inf
    assert reduced_ew(-3.0, 5000.0) == -np.inf
    assert reduced_ew(50.0, -1.0) == -np.inf
    assert reduced_ew(50.0, 0.0) == -np.inf


def test_width_error_budget_composes_in_quadrature():
    """sigma_W^2 = sigma_phot^2 + sigma_cont^2, and sigma_REW = sigma_W / (W ln10).

    Both are identities on values the measurement already holds, so the case asserts the
    composition rather than the size of either term. The two terms are independent by
    construction: the photon term comes from the per-pixel flux noise and the continuum
    term from a fractional offset of the continuum level, and nothing links them.

    Noise model: a constant per-pixel 1 sigma flux error of 1/80 (a signal-to-noise of
    80 per pixel on a unit continuum), and a continuum placement uncertain by 2%.

    The reduced-width relation follows from differentiating REW = log10(W/lambda) at
    fixed lambda: d REW/d W = 1/(W ln 10).

    Tolerance: rel 1e-12, one arithmetic step on stored values.
    """
    from tumult.utils import get_weq

    x, y = gaussian_line()
    sigma = np.full(x.size, 1.0 / 80.0)
    measured = get_weq(x, y, x0=LINE_LAM0, profile="trapz", return_profile=True,
                       sigma=sigma, syst=0.02)

    assert measured.err_phot_maa > 0.0
    assert measured.err_cont_maa > 0.0
    assert measured.width_err_maa == pytest.approx(
        np.hypot(measured.err_phot_maa, measured.err_cont_maa), rel=1e-12)
    assert measured.rew_err == pytest.approx(
        measured.width_err_maa / (measured.maa * np.log(10.0)), rel=1e-12)
    assert measured.rew_err_phot == pytest.approx(
        measured.err_phot_maa / (measured.maa * np.log(10.0)), rel=1e-12)
