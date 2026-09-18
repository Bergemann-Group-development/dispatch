"""Excitation and ionisation balance: recovering the slopes a line set was built with.

Requirement S13 of `reports/test-requirements.md`.

Assumptions:

* per-line abundances are absolute, A(X) in dex, and each carries a 1 sigma error in
  dex. The diagnostics are element agnostic, so the species label only has to name an
  ionisation stage; iron is used because it is what a real balance uses, and nothing
  keys on the symbol;
* excitation potentials are in electron volts and reduced equivalent widths are
  dimensionless. The widths that produce them are in milli-Angstrom and the wavelengths
  in Angstrom on the air scale (the user-facing medium; see *wavelength medium* in
  `docs/GLOSSARY.md`), which is the pairing `reduced_ew` assumes;
* every line carries the same usable 1 sigma error here, so the weighting fallback for a
  line with no usable error is not exercised by the recovery cases and is left to the
  degenerate-input case below;
* the constructed set is in the weak-line regime, where a single abundance per line is
  meaningful. Nothing here synthesises a spectrum: the abundances are the input.
"""
import numpy as np
import pytest

import tumult.sparams as sparams
from tumult.abundance import LineResult

SOLAR_FE = 7.50
LAM0_AA = 5000.0
SIGMA_A = 0.05          # the 1 sigma abundance error every constructed line carries

# Excitation potentials [eV] and reduced widths chosen so that the two abscissae are
# uncorrelated over the set: sum (E - <E>)(R - <R>) = 0. Each univariate slope then
# returns its own coefficient and neither leaks into the other.
E_LOW = np.array([1.0, 2.0, 3.0, 4.0])
REW = np.array([-4.5, -5.5, -5.5, -4.5])
# REW = log10(W/lambda) with W in mAA and lambda in AA, so W = 10^REW * lambda * 1000.
EW_MAA = 10.0**REW * LAM0_AA * 1e3

S_EXC_TRUE = 0.03       # dex per eV
S_REW_TRUE = 0.25       # dex per dex of reduced width
A_BASE = 7.50           # the abundance at the mean excitation and mean reduced width


def neutral_lines(include=(True, True, True, True)):
    """The constructed Fe I set: A = A_BASE + s_exc (E - <E>) + s_rew (R - <R>)."""
    abundances = (A_BASE + S_EXC_TRUE * (E_LOW - E_LOW.mean())
                  + S_REW_TRUE * (REW - REW.mean()))
    return [LineResult(lam0=LAM0_AA, ion="Fe I", e_low=float(E_LOW[i]), loggf=-1.0,
                       ew=float(EW_MAA[i]), a=float(abundances[i]), sigma_a=SIGMA_A,
                       include=bool(include[i]))
            for i in range(len(E_LOW))]


def test_the_excitation_slope_and_its_error_are_recovered():
    """s_exc = 0.03 dex/eV, with sigma = 1/sqrt(S_xx).

    Derivation. The excitation potentials are 1, 2, 3 and 4 eV, so their mean is 2.5 and
    the deviations are -1.5, -0.5, 0.5 and 1.5. Every line has sigma = 0.05, hence a
    weight of 1/0.0025 = 400. The weighted sum of squares is
      S_xx = 400 * (2.25 + 0.25 + 0.25 + 2.25) = 400 * 5 = 2000,
    so the closed-form slope error of a weighted least-squares line is
      sigma_slope = 1/sqrt(2000) = 0.02236068 dex/eV.
    The slope itself is the coefficient the set was built with, because the abundances
    lie exactly on the plane and the reduced-width term is orthogonal to the excitation
    term over this set.

    Tolerance: rel 1e-12. Both are closed forms evaluated on the inputs; the weighting
    cancels in the slope and appears in the error only as the common factor above.
    """
    diagnostics = sparams.diagnostics(neutral_lines(), 0.0, SOLAR_FE)
    assert diagnostics["n_1"] == 4
    assert diagnostics["s_exc"] == pytest.approx(S_EXC_TRUE, rel=1e-12)
    assert diagnostics["s_exc_err"] == pytest.approx(1.0 / np.sqrt(2000.0), rel=1e-12)


def test_the_reduced_width_slope_and_its_error_are_recovered():
    """s_rew = 0.25 dex per dex, with sigma = 1/sqrt(S_xx).

    Derivation. The reduced widths are -4.5, -5.5, -5.5 and -4.5, so their mean is -5
    and the deviations are +-0.5. With the same weight of 400 per line,
      S_xx = 400 * 4 * 0.25 = 400 and sigma_slope = 1/sqrt(400) = 0.05 dex per dex.
    The widths themselves are 158.11388 and 15.811388 mAA at 5000 AA, which is what
    reduced_ew turns back into those reduced widths.

    This is a separate failure mode from the excitation slope: it reaches the same fit
    through the reduced width, so a wrong wavelength or width unit shows up here and
    nowhere else.

    Tolerance: rel 1e-12, closed form as above, with the round trip through log10
    contributing rounding only.
    """
    diagnostics = sparams.diagnostics(neutral_lines(), 0.0, SOLAR_FE)
    assert diagnostics["s_rew"] == pytest.approx(S_REW_TRUE, rel=1e-12)
    assert diagnostics["s_rew_err"] == pytest.approx(1.0 / np.sqrt(400.0), rel=1e-12)


def test_the_ionisation_offset_is_the_difference_between_the_two_stages():
    """d_ion = <A(Fe II)> - <A(Fe I)>, built as an offset of 0.12 dex.

    Derivation. The neutral lines average to A_BASE = 7.50, because both constructed
    trends are taken about their own means and the weights are equal. Three singly
    ionised lines are given 7.62 each, so d_ion = 7.62 - 7.50 = 0.12 dex.

    The metallicity diagnostic is the overall mean against the scaled-solar value:
    the seven lines average to (4*7.50 + 3*7.62)/7 = 7.5514286, and with a solar iron
    abundance of 7.50 and a model [Fe/H] of 0 that leaves
    d_feh = 0.0514286 dex.

    Tolerance: rel 1e-12, weighted means of equal-weight values.
    """
    lines = neutral_lines()
    lines += [LineResult(lam0=5100.0 + i, ion="Fe II", e_low=3.0, loggf=-1.0, ew=40.0,
                         a=A_BASE + 0.12, sigma_a=SIGMA_A, include=True)
              for i in range(3)]

    diagnostics = sparams.diagnostics(lines, 0.0, SOLAR_FE)
    assert diagnostics["n_1"] == 4 and diagnostics["n_2"] == 3
    assert diagnostics["a_1"] == pytest.approx(A_BASE, rel=1e-12)
    assert diagnostics["d_ion"] == pytest.approx(0.12, rel=1e-12)
    assert diagnostics["d_feh"] == pytest.approx(
        (4 * A_BASE + 3 * (A_BASE + 0.12)) / 7 - SOLAR_FE, rel=1e-12)


def test_excluded_and_non_finite_lines_do_not_enter_the_fit():
    """A switched-off line and a line without an abundance leave the slopes untouched.

    Both extra lines are placed far off the constructed trend, so either one entering
    would move the excitation slope by more than a dex per eV. The count is exact,
    because it is a count of what was used.
    """
    lines = neutral_lines()
    lines.append(LineResult(lam0=5200.0, ion="Fe I", e_low=2.0, loggf=-1.0, ew=50.0,
                            a=99.0, sigma_a=SIGMA_A, include=False))
    lines.append(LineResult(lam0=5300.0, ion="Fe I", e_low=2.0, loggf=-1.0, ew=50.0,
                            a=float("nan"), sigma_a=SIGMA_A, include=True))

    diagnostics = sparams.diagnostics(lines, 0.0, SOLAR_FE)
    assert diagnostics["n_1"] == 4
    assert diagnostics["s_exc"] == pytest.approx(S_EXC_TRUE, rel=1e-12)
    assert diagnostics["a_1"] == pytest.approx(A_BASE, rel=1e-12)


def test_a_stage_without_lines_reports_nan_rather_than_a_number():
    """With no singly ionised line there is no ionisation offset to report.

    NaN is the answer a solver reads as "this parameter is undetermined here", so it
    must not be silently replaced by zero, which would read as a satisfied balance.
    """
    diagnostics = sparams.diagnostics(neutral_lines(), 0.0, SOLAR_FE)
    assert diagnostics["n_2"] == 0
    assert np.isnan(diagnostics["d_ion"])


@pytest.mark.parametrize("x, y, sigma, reason", [
    ([1.0, 2.0], [1.0, 2.0], None, "two points cannot support a slope and its error"),
    ([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], None, "a degenerate abscissa has no slope"),
])
def test_a_slope_that_is_not_determined_is_nan(x, y, sigma, reason):
    """Fewer than three points, or one abscissa repeated, gives NaN rather than a fit.

    A weighted least-squares line through two points is exact with no residual and no
    meaningful error, and a vertical set of points has an undefined slope. Returning NaN
    keeps a caller from treating either as a measured trend. Exact, since NaN is a state.
    """
    slope, slope_err, intercept = sparams.weighted_slope(x, y, sigma)
    assert np.isnan(slope), reason
    assert np.isnan(slope_err)
    assert np.isnan(intercept)


def test_a_line_without_a_usable_error_falls_back_instead_of_poisoning_the_fit():
    """A missing or non-positive sigma takes the median of the usable ones.

    Derivation: the four points lie exactly on y = x, so any set of positive weights
    returns a slope of 1. What is being checked is that a NaN or a zero error does not
    propagate into the weighted sums, which would make the whole fit NaN and silently
    remove the diagnostic from a solve.

    Tolerance: rel 1e-12, since the data are exactly collinear.
    """
    slope, _err, intercept = sparams.weighted_slope(
        [0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0], [0.1, float("nan"), 0.1, 0.0])
    assert slope == pytest.approx(1.0, rel=1e-12)
    assert intercept == pytest.approx(0.0, abs=1e-12)
