"""Guards the Phase 1 utils/auxfuncs consolidation contract.

Originally (pre-Phase 1) this file proved the duplicated implementations in
utils.py and auxfuncs.py were numerically identical; after the merge it asserts the
stronger property that the auxfuncs names ARE the utils/physics objects (aliases),
plus the deliberate divergences that were kept apart.
"""
import numpy as np
import pytest

from tumult import auxfuncs, physics, utils


RNG = np.random.default_rng(20260703)
XX = np.sort(RNG.uniform(4000.0, 4010.0, 50))
YY = 1 - 0.5 * np.exp(-0.5 * ((XX - 4005.0) / 0.8) ** 2) + 0.01 * RNG.standard_normal(50)


def test_single_implementations():
    """auxfuncs re-exports the canonical utils/physics implementations - same objects."""
    for name in ("pchip_mono", "parabolic_min", "LinArrayInterp",
                 "linear_edge_normalisation", "get_bisector"):
        assert getattr(auxfuncs, name) is getattr(utils, name), name
    for name in ("vac2air", "air2vac", "planck", "voigt_function"):
        assert getattr(auxfuncs, name) is getattr(physics, name), name


def test_constants_single_source():
    from tumult import const

    for mod in (auxfuncs,):
        for name in ("EE", "HH", "CC", "BK", "EM", "UU", "PI", "HCE", "HC2", "HCK",
                     "EK", "XMH", "SBOL", "CC_KMS"):
            assert getattr(mod, name) == getattr(const, name), name
    import tumult.atom as tatom
    assert tatom.saha_const is const.saha_const


def test_lowpass_implementations_DIVERGE():
    """utils.lowpass (Gaussian roll-off at resolution R in ln-lambda) and the legacy
    brick-wall FFT cut are different filters; the legacy one is now
    auxfuncs.lowpass_fft_thresh, with `auxfuncs.lowpass` as a deprecated alias."""
    a = utils.lowpass(XX, YY, 5000.0)
    b = auxfuncs.lowpass_fft_thresh(XX, YY, 15)
    assert a.shape == b.shape == YY.shape
    assert not np.allclose(a, b)
    with pytest.warns(DeprecationWarning):
        c = auxfuncs.lowpass(XX, YY, 15)
    np.testing.assert_array_equal(b, c)


def test_voigt_profile_implementations_DIVERGE():
    """utils.voigt_profile(x, x0, amp, sigma, rho) delegates to VoigtProfile.model
    (peak-normalised); auxfuncs.voigt_profile(xx, x0, a, sigma, amp) is the legacy
    un-normalised parameterisation. They are NOT the same model and are kept apart."""
    x = np.linspace(4999, 5001, 201)
    modern = utils.voigt_profile(x, 5000.0, 0.4, 0.1, 0.2)
    legacy = auxfuncs.voigt_profile(x, 5000.0, 0.2, 0.1, 0.4)
    assert not np.allclose(modern, legacy)


def test_calc_teff_lives_in_utils_only():
    """The raw-CGS auxfuncs.calc_teff (dead planck_thresh branch) was deleted; the
    unit-aware utils version is the single implementation and must recover a
    blackbody temperature."""
    import astropy.units as u
    from astropy.modeling.models import BlackBody

    assert not hasattr(auxfuncs, "calc_teff")
    teff_in = 5500.0
    lam = np.geomspace(100, 5e5, 20000) * u.AA
    bb = BlackBody(teff_in * u.K, scale=1 * u.erg / (u.cm**2 * u.AA * u.s * u.sr))
    flux = np.pi * bb(lam) * u.sr
    assert utils.calc_teff(lam, flux).to(u.K).value == pytest.approx(teff_in, rel=1e-3)


def test_pad_boolean_defined_once():
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(auxfuncs))
    defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "pad_boolean"]
    assert len(defs) == 1
