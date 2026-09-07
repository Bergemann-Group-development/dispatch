"""Pin the equivalent-width machinery: every registered profile, all return conventions.

The EW path is physics-critical (abundances hang off it). These tests freeze the current
behaviour of ``Profile.measure`` / ``utils.get_weq`` / ``Spectrum.get_weq`` so the Phase 3
consolidation onto a single Profile-based path can prove it changed nothing.
"""
import astropy.units as u
import numpy as np
import pytest

from cases import LINE_AMP, LINE_SIG, LINE_X0, ew_profile_cases, gauss_line
from conftest import assert_matches_golden, load_golden


# The continuum-placement error columns of the two FITTED profiles come out of the
# least-squares covariance, whose last digits track the SciPy/BLAS build rather than
# anything here: they drift ~1e-9 across environments while the widths themselves stay
# exact. The EWs, the fit parameters and every integrated profile keep the 1e-12 default.
_FIT_COV_RTOL = {"*_err_cont_mAA": 1e-8, "*_rew_err": 1e-8}


def test_matches_golden():
    assert_matches_golden(ew_profile_cases(), load_golden("ew_profiles"),
                          rtol_overrides=_FIT_COV_RTOL)


def test_registry_is_complete():
    from tumult.profile import profile_names

    assert set(profile_names()) == {"gauss", "voigt", "asym_gauss", "trapz", "pchip"}


def test_gauss_ew_matches_analytic():
    """W = amp*sqrt(2*pi)*sig for a Gaussian - the fitted and integrated EWs must agree."""
    from tumult.utils import get_weq

    x, y = gauss_line()
    analytic_maa = LINE_AMP * np.sqrt(2 * np.pi) * LINE_SIG * 1e3
    for profile in ("gauss", "trapz", "pchip", "voigt", "asym_gauss"):
        prof = get_weq(x, y, x0=LINE_X0, profile=profile, return_profile=True)
        assert prof.mAA == pytest.approx(analytic_maa, rel=1e-3), profile


def test_spectrum_get_weq_equals_utils_get_weq():
    """Spectrum.get_weq is a passthrough; the two entry points must agree exactly."""
    from tumult.spectrum import Spectrum
    from tumult.utils import get_weq

    x, y = gauss_line()
    sp = Spectrum(x, y)
    for profile in ("gauss", "trapz"):
        direct = get_weq(x, y, x0=LINE_X0, profile=profile, return_profile=True)
        via_spec = sp.get_weq(x0=LINE_X0 * u.AA, profile=profile, return_profile=True)
        assert via_spec.mAA == pytest.approx(direct.mAA, rel=1e-12), profile
        assert via_spec.rew == pytest.approx(direct.rew, rel=1e-12), profile


def test_error_budget_composition():
    """Photon and continuum width errors combine in quadrature to the total."""
    from tumult.utils import get_weq

    x, y = gauss_line()
    sigma = np.full(x.size, 1.0 / 80.0)
    prof = get_weq(x, y, x0=LINE_X0, profile="trapz", return_profile=True,
                   sigma=sigma, syst=0.02)
    assert prof.err_phot_mAA > 0
    assert prof.err_cont_mAA > 0
    assert prof.width_err_mAA == pytest.approx(
        np.hypot(prof.err_phot_mAA, prof.err_cont_mAA), rel=1e-12)
    # sigma_REW = sigma_W / (W ln10)
    assert prof.rew_err == pytest.approx(
        prof.width_err_mAA / (prof.mAA * np.log(10)), rel=1e-12)


def test_reduced_ew_guards():
    from tumult.profile import reduced_ew

    assert reduced_ew(50.0, 5000.0) == pytest.approx(np.log10(50e-3 / 5000.0))
    assert reduced_ew(0.0, 5000.0) == -np.inf
    assert reduced_ew(-3.0, 5000.0) == -np.inf
    assert reduced_ew(50.0, -1.0) == -np.inf
