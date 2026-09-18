"""Abundances: the element mean over lines, and the inversion of a curve of growth.

Requirements S4 and S5 of `reports/test-requirements.md`.

Shared assumptions:

* abundances are absolute, A(X) = log10(N_X/N_H) + 12, in dex; an error on one is a
  1 sigma error in dex;
* equivalent widths are in milli-Angstrom and wavelengths in Angstrom on the air scale
  (the user-facing medium; see *wavelength medium* in `docs/GLOSSARY.md`). The curve of
  growth is a relation between two numbers, so the medium only has to be consistent
  between the width and the wavelength it is reduced by;
* the lines are measured on a continuum-normalised spectrum, which is what makes an
  equivalent width the area of the profile below a continuum of exactly 1.0;
* the weak-line regime: the constructed curve of growth is log-linear in
  (REW, A), which is the unsaturated part of a real curve. Saturation is out of scope
  here and would break the closed forms below.
"""
from types import SimpleNamespace

import numpy as np
import pytest

# The constructed line set of S4. Each line carries an abundance, a total 1 sigma error
# and a photon/continuum split obeying sigma^2 = sigma_phot^2 + sigma_cont^2 exactly
# (3-4-5 triples scaled), so the quadrature decomposition is not an approximation here.
LINES = [
    SimpleNamespace(a=7.50, sigma_a=0.10, sigma_a_phot=0.06, sigma_a_cont=0.08,
                    include=True),
    SimpleNamespace(a=7.40, sigma_a=0.20, sigma_a_phot=0.12, sigma_a_cont=0.16,
                    include=True),
    SimpleNamespace(a=7.60, sigma_a=0.20, sigma_a_phot=0.16, sigma_a_cont=0.12,
                    include=True),
    # switched off by the user: its wild abundance must not move the mean at all
    SimpleNamespace(a=9.00, sigma_a=0.10, sigma_a_phot=0.06, sigma_a_cont=0.08,
                    include=False),
    # carries no usable error, so it has zero weight and cannot enter an
    # inverse-variance mean
    SimpleNamespace(a=6.00, sigma_a=np.inf, sigma_a_phot=np.inf, sigma_a_cont=np.inf,
                    include=True),
]


def test_element_mean_is_the_inverse_variance_weighted_mean():
    """Hand-derived mean and formal error over the three included, weighted lines.

    Derivation. Weights w = 1/sigma^2 are 100, 25 and 25, summing to 150.
    Mean = (100*7.50 + 25*7.40 + 25*7.60)/150 = (750 + 185 + 190)/150 = 1125/150 = 7.50.
    Formal error = sqrt(1/150) = 0.08164966.
    Weighted population scatter = sqrt(sum w (a - mean)^2 / sum w)
      = sqrt((100*0 + 25*0.01 + 25*0.01)/150) = sqrt(0.5/150) = sqrt(1/300) = 0.05773503.
    The scatter is the population form, with no sample (n-1) correction: it is a
    cross-check on the formal error, deliberately conservative.

    Tolerance: rel 1e-12. Every number above is closed-form arithmetic on the inputs,
    and the count is exact because it is a count.
    """
    from tumult.abundance import aggregate

    result = aggregate(LINES)
    assert result["n"] == 3
    assert result["a_mean"] == pytest.approx(7.50, rel=1e-12)
    assert result["sigma_mean"] == pytest.approx(np.sqrt(1.0 / 150.0), rel=1e-12)
    assert result["sigma_scatter"] == pytest.approx(np.sqrt(1.0 / 300.0), rel=1e-12)


def test_photon_and_continuum_shares_add_to_the_formal_error():
    """The mean's error splits into the two shares that build each line's error.

    Derivation. With weights w_i the mean is sum w_i a_i / sum w_i, so a component k
    contributes sum w_i^2 sigma_{k,i}^2 / (sum w_i)^2 to its variance. With
    sum w = 150 and denominator 22500:
      photon: (100^2*0.06^2 + 25^2*0.12^2 + 25^2*0.16^2)/22500
              = (36 + 9 + 16)/22500 = 61/22500 = 0.002711111, so 0.05206833;
      continuum: (100^2*0.08^2 + 25^2*0.16^2 + 25^2*0.12^2)/22500
              = (64 + 16 + 9)/22500 = 89/22500 = 0.003955556, so 0.06289320.
    Their quadrature sum is (61 + 89)/22500 = 1/150, the formal variance above. That
    identity holds because each constructed line satisfies
    sigma^2 = sigma_phot^2 + sigma_cont^2 exactly; for measured lines with a different
    split the decomposition would only be approximate.

    Tolerance: rel 1e-12, closed-form arithmetic throughout.
    """
    from tumult.abundance import aggregate

    result = aggregate(LINES)
    assert result["sigma_phot"] == pytest.approx(np.sqrt(61.0 / 22500.0), rel=1e-12)
    assert result["sigma_cont"] == pytest.approx(np.sqrt(89.0 / 22500.0), rel=1e-12)
    assert np.hypot(result["sigma_phot"], result["sigma_cont"]) == pytest.approx(
        result["sigma_mean"], rel=1e-12)


def test_excluded_and_zero_weight_lines_do_not_enter_the_mean():
    """Dropping the two lines that must not count leaves the same answer.

    An excluded line is a user decision and a line without a usable error has no place
    in an inverse-variance mean; either one entering would move the mean, since both sit
    far from it (9.00 and 6.00 against 7.50). Comparing against the mean of the three
    weighted lines alone is exact: the same arithmetic on the same inputs.
    """
    from tumult.abundance import aggregate

    reduced = aggregate(LINES[:3])
    full = aggregate(LINES)
    assert full["n"] == reduced["n"] == 3
    assert full["a_mean"] == reduced["a_mean"]
    assert full["sigma_mean"] == reduced["sigma_mean"]


def test_an_empty_line_set_gives_nan_with_zero_lines():
    """No line is not an abundance of zero: the mean is NaN and the count is 0.

    Exact on the count, and NaN asserted as NaN, because both are discrete outcomes of
    a decision rather than measured quantities. A caller distinguishes "not measured"
    from "measured as zero" on this.
    """
    from tumult.abundance import aggregate

    empty = aggregate([])
    assert empty["n"] == 0
    assert np.isnan(empty["a_mean"])
    assert np.isnan(empty["sigma_mean"])
    assert np.isnan(empty["sigma_scatter"])

    all_excluded = aggregate([SimpleNamespace(a=7.5, sigma_a=0.1, sigma_a_phot=0.06,
                                              sigma_a_cont=0.08, include=False)])
    assert all_excluded["n"] == 0
    assert np.isnan(all_excluded["a_mean"])


# --- curve of growth ------------------------------------------------------------------

# A curve exactly log-linear in (REW, A): W(A) = W0 * 10^(k A) with W0 = 55 mAA and
# k = 0.62 dex of width per dex of abundance, on seven rungs spanning +-0.6 dex.
COG_W0_MAA = 55.0
COG_SLOPE = 0.62
COG_ABUNDANCES = np.linspace(-0.6, 0.6, 7)
COG_WIDTHS_MAA = COG_W0_MAA * 10.0 ** (COG_SLOPE * COG_ABUNDANCES)
LAM0_AA = 5000.0


def test_inverting_a_log_linear_curve_of_growth_returns_the_analytic_abundance():
    """A = log10(W/W0)/k on a curve that is a straight line in (REW, A).

    Derivation. W(A) = 55 * 10^(0.62 A) mAA inverts to A = log10(W/55)/0.62. For a
    measured 50 mAA that is log10(50/55)/0.62 = -0.04139269/0.62 = -0.06676238 dex.
    The local slope is dA/dREW = 1/k = 1/0.62 = 1.61290323 dex per dex, constant over
    the whole curve because the curve is straight.
    The measured error propagates as sigma_REW = sigma_W/(W ln10), so with
    sigma_W = 2.5 mAA on 50 mAA, sigma_REW = 2.5/(50 * 2.30258509) = 0.02171472 and
    sigma_A = 0.02171472 * 1.61290323 = 0.03502375 dex.
    The reduced width itself is log10(0.050/5000) = -5 exactly, since lam0 is supplied
    in Angstrom and the widths in milli-Angstrom.

    Tolerance: rel 1e-10. The shape-preserving interpolant through collinear rungs is
    that line, so the only error is rounding through two logarithms and the spline's
    coefficient arithmetic. Anything larger would mean the interpolant is bending
    where the data do not.
    """
    from tumult.abundance import cog_invert

    inverted = cog_invert(COG_ABUNDANCES, COG_WIDTHS_MAA, ew_meas=50.0, ew_err=2.5,
                          lam0_aa=LAM0_AA)
    assert inverted["a"] == pytest.approx(np.log10(50.0 / COG_W0_MAA) / COG_SLOPE,
                                          rel=1e-10)
    assert inverted["dadrew"] == pytest.approx(1.0 / COG_SLOPE, rel=1e-10)
    assert inverted["rew"] == pytest.approx(-5.0, rel=1e-10)
    assert inverted["sigma_a"] == pytest.approx(
        2.5 / (50.0 * np.log(10.0)) / COG_SLOPE, rel=1e-10)


def test_a_measurement_inside_the_synthesised_span_is_not_flagged():
    """50 mAA lies between the weakest and strongest rung, so the answer is interpolated.

    The rungs span 55 * 10^(+-0.372) = 23.35 to 129.53 mAA. The flag is a discrete
    decision and is asserted exactly.
    """
    from tumult.abundance import cog_invert

    inverted = cog_invert(COG_ABUNDANCES, COG_WIDTHS_MAA, ew_meas=50.0,
                          lam0_aa=LAM0_AA)
    assert inverted["out_of_range"] is False


@pytest.mark.parametrize("ew_meas", [20.0, 150.0])
def test_a_measurement_outside_the_synthesised_span_is_flagged(ew_meas):
    """Below 23.35 mAA and above 129.53 mAA the answer is an extrapolation.

    Both sides are exercised because they are different failure modes: a line too weak
    for the ladder and a line too strong for it reach the flag through opposite
    comparisons, and an inverted comparison would pass with only one of them.
    """
    from tumult.abundance import cog_invert

    inverted = cog_invert(COG_ABUNDANCES, COG_WIDTHS_MAA, ew_meas=ew_meas,
                          lam0_aa=LAM0_AA)
    assert inverted["out_of_range"] is True
