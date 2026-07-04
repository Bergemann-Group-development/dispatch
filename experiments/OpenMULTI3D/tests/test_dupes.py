"""TEMPORARY (delete at the end of Phase 1): documents the utils.py / auxfuncs.py
duplications so the merge is provably behaviour-preserving.

Package-level note: tumult/__init__.py wildcard-imports auxfuncs BEFORE utils, so for
every duplicated name the utils version already wins at package level - consolidating on
utils preserves `tumult.<name>` behaviour by construction.
"""
import numpy as np
import pytest

from tumult import auxfuncs, utils


RNG = np.random.default_rng(20260703)
XX = np.sort(RNG.uniform(4000.0, 4010.0, 50))
YY = 1 - 0.5 * np.exp(-0.5 * ((XX - 4005.0) / 0.8) ** 2) + 0.01 * RNG.standard_normal(50)


def test_pchip_mono_identical():
    grid = np.linspace(XX[0], XX[-1], 500)
    np.testing.assert_array_equal(utils.pchip_mono(XX, YY)(grid),
                                  auxfuncs.pchip_mono(XX, YY)(grid))


def test_parabolic_min_identical():
    assert utils.parabolic_min(XX, YY) == auxfuncs.parabolic_min(XX, YY)


def test_linarrayinterp_identical():
    newx = np.linspace(XX[2], XX[-3], 77)
    np.testing.assert_array_equal(utils.LinArrayInterp(XX, YY, newx),
                                  auxfuncs.LinArrayInterp(XX, YY, newx))


def test_linear_edge_normalisation_identical():
    np.testing.assert_array_equal(utils.linear_edge_normalisation(XX, YY),
                                  auxfuncs.linear_edge_normalisation(XX, YY))


def test_get_bisector_identical():
    ux, uy = utils.get_bisector(XX, YY, x0=4005.0)
    ax, ay = auxfuncs.get_bisector(XX, YY, x0=4005.0)
    np.testing.assert_array_equal(ux, ax)
    np.testing.assert_array_equal(uy, ay)


def test_voigt_function_identical():
    a = np.array([0.0, 0.1, 0.5])
    v = np.linspace(-4, 4, 41)
    for ai in a:
        np.testing.assert_array_equal(utils.voigt_function(ai, v),
                                      auxfuncs.voigt_function(ai, v))


def test_lowpass_implementations_DIVERGE():
    """The two `lowpass` functions are different filters with different signatures:
    utils.lowpass(xx, yy, resolution) - Gaussian roll-off in ln-lambda space;
    auxfuncs.lowpass(xx, yy, thresh) - brick-wall FFT cut on a linear grid.
    Nothing in the codebase calls the auxfuncs one; Phase 1 renames it to
    lowpass_fft_thresh instead of merging."""
    a = utils.lowpass(XX, YY, 5000.0)
    b = auxfuncs.lowpass(XX, YY, 5000.0)
    assert a.shape == b.shape == YY.shape
    assert not np.allclose(a, b)


def test_voigt_profile_implementations_DIVERGE():
    """utils.voigt_profile(x, x0, amp, sigma, rho) delegates to VoigtProfile.model
    (peak-normalised); auxfuncs.voigt_profile(xx, x0, a, sigma, amp) is an older
    un-normalised parameterisation with a different argument order. They are NOT the
    same model - Phase 1 keeps them apart (profile.py is authoritative)."""
    x = np.linspace(4999, 5001, 201)
    modern = utils.voigt_profile(x, 5000.0, 0.4, 0.1, 0.2)
    legacy = auxfuncs.voigt_profile(x, 5000.0, 0.2, 0.1, 0.4)
    assert not np.allclose(modern, legacy)


def test_calc_teff_conventions():
    """utils.calc_teff works in astropy units, auxfuncs.calc_teff in raw CGS (and its
    planck_thresh branch is dead - it calls an undefined `crop_spectrum`). Both must
    recover the temperature of a blackbody in their own convention."""
    import astropy.units as u
    from astropy.modeling.models import BlackBody

    teff_in = 5500.0
    lam = np.geomspace(100, 5e5, 20000) * u.AA
    bb = BlackBody(teff_in * u.K, scale=1 * u.erg / (u.cm**2 * u.AA * u.s * u.sr))
    flux = np.pi * bb(lam) * u.sr
    teff_u = utils.calc_teff(lam, flux)
    assert teff_u.to(u.K).value == pytest.approx(teff_in, rel=1e-3)

    lam_cm = lam.to(u.cm).value
    flux_cgs = flux.to(u.erg / (u.cm**2 * u.cm * u.s)).value
    teff_a = auxfuncs.calc_teff(lam_cm, flux_cgs)
    assert teff_a == pytest.approx(teff_in, rel=1e-3)

    with pytest.raises(NameError):  # the dead branch, documented
        auxfuncs.calc_teff(lam_cm, flux_cgs, planck_thresh=1e-4)


def test_pad_boolean_defined_twice_identically():
    """auxfuncs defines pad_boolean twice (~695 and ~848); the second silently shadows
    the first. They are AST-identical today - Phase 1 deletes one."""
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(auxfuncs))
    defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "pad_boolean"]
    assert len(defs) == 2, "Phase 1 done? update/delete this test"
    assert ast.dump(defs[0]) == ast.dump(defs[1])
