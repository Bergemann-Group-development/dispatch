"""Pin the excitation–ionization balance numerics: the weighted slope, the four balance
diagnostics, the damped quasi-Newton solver (with secant Jacobian refinement, step clamps,
grid bounds and held parameters), the uncertainty propagation, and the measurement object's
session round-trip.

The numerics are element-general — the tests use Fe, whose sensitivities seed the solver for
every element, and one Ti case to pin that nothing is keyed to the symbol.
"""
import numpy as np
import pytest

from tumult import measurement as ms
from tumult import sparams as sp
from tumult.abundance import LineResult
from tumult.abundance import reduced_ew

SOLAR_FE = 7.50


def _line(lam0, ion, e_low, ew, a, sigma_a=0.05, include=True):
    return LineResult(lam0=lam0, ion=ion, e_low=e_low, loggf=-1.0, ew=ew,
                      a=a, sigma_a=sigma_a, include=include)


# --- weighted_slope -----------------------------------------------------------------------

def test_weighted_slope_recovers_line():
    x = np.linspace(0.0, 5.0, 20)
    y = 1.3 * x - 0.7
    slope, err, intercept = sp.weighted_slope(x, y, np.full_like(x, 0.1))
    assert slope == pytest.approx(1.3, abs=1e-12)
    assert intercept == pytest.approx(-0.7, abs=1e-12)
    assert err > 0


def test_weighted_slope_weights_matter():
    """A discrepant point with a huge sigma must barely move the fit."""
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = 2.0 * x
    y_off = y.copy()
    y_off[-1] += 10.0
    s = np.array([0.1, 0.1, 0.1, 0.1, 100.0])
    slope, _, _ = sp.weighted_slope(x, y_off, s)
    assert slope == pytest.approx(2.0, abs=0.01)


def test_weighted_slope_degenerate():
    nan = float("nan")
    assert np.isnan(sp.weighted_slope([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])[0])
    assert np.isnan(sp.weighted_slope([1.0, 2.0], [1.0, 2.0])[0])
    # non-finite sigmas fall back rather than poisoning the fit
    slope, _, _ = sp.weighted_slope([0, 1, 2, 3], [0, 1, 2, 3], [0.1, nan, 0.1, 0.0])
    assert slope == pytest.approx(1.0, abs=1e-9)


# --- diagnostics ----------------------------------------------------------------------------

def test_diagnostics_recovers_construction():
    """Fe I lines built with a known E_low slope; Fe II offset by a known ΔA."""
    feh = -0.5
    a0 = SOLAR_FE + feh          # consistency-perfect mean at zero slope
    s_true = 0.04                # dex/eV
    d_true = 0.12
    lines = [_line(5000 + i, "Fe I", e, ew=50.0, a=a0 + s_true * (e - 2.5))
             for i, e in enumerate(np.linspace(0.5, 4.5, 12))]
    lines += [_line(6000 + i, "Fe II", 3.0, ew=40.0, a=a0 + d_true) for i in range(4)]
    d = sp.diagnostics(lines, feh, SOLAR_FE)
    assert d["s_exc"] == pytest.approx(s_true, abs=1e-9)
    assert d["n_1"] == 12 and d["n_2"] == 4
    assert d["a_1"] == pytest.approx(a0, abs=1e-9)        # slope pivots about the mean
    assert d["d_ion"] == pytest.approx(d_true, abs=1e-9)
    # overall mean sits between the stages; d_feh measures it against solar+feh
    assert d["d_feh"] == pytest.approx(d["a_mean"] - (SOLAR_FE + feh), abs=1e-12)


def test_diagnostics_rew_slope():
    """A(Fe I) rising with reduced EW must show as a positive s_rew."""
    ews = np.linspace(20.0, 120.0, 10)
    rews = [reduced_ew(w, 5500.0) for w in ews]
    lines = [_line(5500.0, "Fe I", 2.0, ew=w, a=7.5 + 0.3 * (r - np.mean(rews)))
             for w, r in zip(ews, rews)]
    d = sp.diagnostics(lines, 0.0, SOLAR_FE)
    assert d["s_rew"] == pytest.approx(0.3, abs=1e-6)


def test_diagnostics_excluded_and_nan_lines_ignored():
    lines = [_line(5000 + i, "Fe I", e, 50.0, 7.5) for i, e in enumerate([1.0, 2.0, 3.0, 4.0])]
    lines.append(_line(5100.0, "Fe I", 2.0, 50.0, 99.0, include=False))
    lines.append(_line(5200.0, "Fe I", 2.0, 50.0, float("nan")))
    d = sp.diagnostics(lines, 0.0, SOLAR_FE)
    assert d["n_1"] == 4
    assert d["a_1"] == pytest.approx(7.5, abs=1e-12)
    assert np.isnan(d["d_ion"])          # no Fe II lines → NaN, solver holds logg


def test_diagnostics_are_element_agnostic():
    """Nothing is keyed to Fe: the same construction in Ti reads back identically."""
    solar_ti = 4.97
    lines = [_line(5000 + i, "Ti I", e, 50.0, solar_ti - 1.0 + 0.04 * (e - 2.5))
             for i, e in enumerate(np.linspace(0.5, 4.5, 8))]
    lines += [_line(6000.0, "Ti II", 3.0, 40.0, solar_ti - 1.0 + 0.1)]
    d = sp.diagnostics(lines, -1.0, solar_ti)
    assert d["n_1"] == 8 and d["n_2"] == 1
    assert d["s_exc"] == pytest.approx(0.04, abs=1e-9)
    assert d["d_ion"] == pytest.approx(0.1, abs=1e-9)
    assert d["d_feh"] == pytest.approx(d["a_mean"] - (solar_ti - 1.0), abs=1e-12)


def test_converged():
    good = {"s_exc": 0.001, "s_rew": 0.005, "d_ion": 0.005, "d_feh": -0.005}
    assert sp.converged(good)
    assert not sp.converged({**good, "d_ion": 0.05})
    # a NaN diagnostic doesn't block convergence of the others…
    assert sp.converged({**good, "d_ion": float("nan")})
    # …but an all-NaN set is not converged
    assert not sp.converged({k: float("nan") for k in sp.DIAG_KEYS})


def test_converged_ignores_held_parameters():
    """A held parameter's diagnostic is reported but cannot gate: it can never be nulled."""
    off = {"s_exc": 0.001, "s_rew": 0.005, "d_ion": 0.5, "d_feh": -0.005}
    assert not sp.converged(off)                      # the ionisation offset is way out…
    free = {"teff": True, "logg": False, "feh": True, "vmic": True}
    assert sp.converged(off, free=free)               # …but log g was held there
    assert set(sp.free_diags(free)) == {"s_exc", "s_rew", "d_feh"}
    # holding everything is a diagnostic evaluation, not a solution
    assert not sp.converged(off, free={k: False for k in sp.PARAM_KEYS})


# --- solver ---------------------------------------------------------------------------------

TRUTH = {"teff": 5900.0, "logg": 4.30, "feh": -0.20, "vmic": 1.20}
# a fabricated diagonal linear model diag = J_true * (p - p_true), with sensitivities
# deliberately different from sparams.DEFAULT_JAC so the secant refinement has work to do
J_TRUE = {"teff": -8.0e-5, "logg": 0.55, "feh": -1.0, "vmic": -0.22}


def _model_diags(params):
    d = {}
    for key in sp.PARAM_KEYS:
        d[sp.DIAG_FOR[key]] = J_TRUE[key] * (params[key] - TRUTH[key])
    return d


def test_step_converges_on_linear_model():
    params = {"teff": 5400.0, "logg": 4.8, "feh": 0.2, "vmic": 0.8}
    state = sp.SolverState(params=dict(params))
    for _ in range(15):
        diags = _model_diags(state.params)
        if sp.converged(diags):
            break
        sp.step(state, diags)
    assert sp.converged(_model_diags(state.params))
    assert state.params["teff"] == pytest.approx(TRUTH["teff"], abs=25.0)
    assert state.params["logg"] == pytest.approx(TRUTH["logg"], abs=0.03)
    assert state.params["feh"] == pytest.approx(TRUTH["feh"], abs=0.02)
    assert state.params["vmic"] == pytest.approx(TRUTH["vmic"], abs=0.05)


def test_step_respects_clamps_and_bounds():
    params = {"teff": 4000.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0}
    state = sp.SolverState(params=dict(params),
                           bounds={"teff": (3900.0, 4150.0), "feh": (-0.05, 0.05)})
    diags = {"s_exc": 0.5, "s_rew": 0.0, "d_ion": 0.0, "d_feh": -3.0}
    new = sp.step(state, diags)
    # a huge excitation slope is clamped to one STEP_CLAMP, then to the grid box
    assert new["teff"] <= 4150.0
    assert new["teff"] - params["teff"] <= sp.STEP_CLAMP["teff"] + 1e-9
    assert -0.05 <= new["feh"] <= 0.05
    # untouched axes hold still on zero diagnostics
    assert new["logg"] == params["logg"] and new["vmic"] == params["vmic"]


def test_step_holds_fixed_parameters():
    """A held parameter does not move, however large its diagnostic."""
    params = {"teff": 5400.0, "logg": 4.8, "feh": 0.2, "vmic": 0.8}
    state = sp.SolverState(params=dict(params))
    free = {"teff": True, "logg": False, "feh": False, "vmic": True}
    new = sp.step(state, {"s_exc": 0.05, "s_rew": 0.1, "d_ion": 0.4, "d_feh": 0.4},
                  free=free)
    assert new["logg"] == params["logg"] and new["feh"] == params["feh"]
    assert new["teff"] != params["teff"] and new["vmic"] != params["vmic"]


def test_step_holds_param_on_nan_diagnostic():
    state = sp.SolverState(params={"teff": 5777.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0})
    new = sp.step(state, {"s_exc": float("nan"), "s_rew": 0.0,
                          "d_ion": 0.1, "d_feh": 0.0})
    assert new["teff"] == 5777.0
    assert new["logg"] != 4.4


def test_secant_updates_keep_sign_and_window():
    state = sp.SolverState(params={"teff": 5400.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0})
    sp.step(state, _model_diags(state.params))
    sp.step(state, _model_diags(state.params))
    lo, hi = sp.JAC_CLAMP
    for key in sp.PARAM_KEYS:
        ref = sp.DEFAULT_JAC[key]
        assert state.jac[key] * ref > 0                      # sign preserved
        assert lo * abs(ref) - 1e-12 <= abs(state.jac[key]) <= hi * abs(ref) + 1e-12
    # the teff entry actually moved toward the true (stronger) sensitivity
    assert abs(state.jac["teff"]) > abs(sp.DEFAULT_JAC["teff"])


def test_uncertainties():
    diags = {"s_exc_err": 0.004, "s_rew_err": 0.02, "d_ion_err": 0.05,
             "d_feh_err": 0.03}
    sig = sp.uncertainties(diags, dict(J_TRUE))
    assert sig["teff"] == pytest.approx(0.004 / 8.0e-5)
    assert sig["logg"] == pytest.approx(0.05 / 0.55)
    assert sig["feh"] == pytest.approx(0.03)
    assert sig["vmic"] == pytest.approx(0.02 / 0.22)
    assert np.isnan(sp.uncertainties({}, dict(J_TRUE))["teff"])
    # a held parameter was not determined here, so it gets no error of its own
    held = sp.uncertainties(diags, dict(J_TRUE),
                            {"teff": True, "logg": False, "feh": True, "vmic": True})
    assert np.isnan(held["logg"]) and held["teff"] == pytest.approx(0.004 / 8.0e-5)


# --- the measurement object -------------------------------------------------------------

def _measurement(**kw):
    lines = [_line(5000.0, "Fe I", 1.0, 50.0, 7.4), _line(5100.0, "Fe II", 3.0, 30.0, 7.6)]
    base = dict(lines=lines, element="Fe", maskset_uid="3", maskset_name="fe.txt",
                maskset_key="4999.6000-5000.4000", ref_name="HD140283",
                params=dict(TRUTH), params0={"teff": 5777.0, "logg": 4.4,
                                             "feh": 0.0, "vmic": 1.0},
                entry="gdwarf", free={k: True for k in sp.PARAM_KEYS},
                history=[{"params": dict(TRUTH), "diags": {"s_exc": 0.0}}],
                sigmas={"teff": 40.0}, diags={"s_exc": 0.0}, jac=dict(J_TRUE),
                converged=True, solar_x=SOLAR_FE, stage=ms.STAGE_SOLVED,
                profile="gauss", snr=100.0, cont_pct=0.5, uid="abc123")
    base.update(kw)
    return ms.Measurement(**base)


def test_measurement_roundtrip():
    res = _measurement(name="Fe bal")
    back = ms.Measurement.from_dict(res.to_dict())
    assert back.params == res.params and back.jac == res.jac
    assert back.converged and back.stage == ms.STAGE_SOLVED
    assert back.solar_x == SOLAR_FE
    assert (back.element, back.ref_name, back.entry) == ("Fe", "HD140283", "gdwarf")
    assert back.free == res.free and back.maskset_key == res.maskset_key
    assert len(back.lines) == 2
    assert back.lines[0].ion == "Fe I"
    assert back.lines[0].a == pytest.approx(7.4)
    assert back.history[0]["diags"]["s_exc"] == 0.0


def test_measurement_counts_and_labels():
    res = _measurement()
    assert res.counts() == (1, 1) and res.has_both_stages()
    assert res.label().startswith("Teff 5900")
    res.converged = False
    assert "not converged" in res.label()
    # every parameter held: a diagnostic evaluation, and the label says so
    held = _measurement(free={k: False for k in sp.PARAM_KEYS})
    assert held.diagnose_only and "diagnostics" in held.label()
    # before a fit, the census is the name
    fresh = _measurement(stage=ms.STAGE_MEASURED, params={})
    assert fresh.label() == "Fe balance · 1 Fe I / 1 Fe II · HD140283"
    fresh.lines[1].include = False
    assert fresh.counts() == (1, 0) and not fresh.has_both_stages()


def test_measurement_reaggregates_without_synthesis():
    """Toggling a line off re-reads the diagnostics from the stored per-line abundances."""
    res = _measurement()
    res.reaggregate()
    assert res.diags["n_1"] == 1 and res.diags["n_2"] == 1
    res.lines[1].include = False
    res.reaggregate()
    assert res.diags["n_2"] == 0 and np.isnan(res.diags["d_ion"])
    assert res.diags["a_1"] == pytest.approx(7.4)


def test_measurement_staleness_follows_the_windows():
    class _Mask:
        def __init__(self, lmin, lmax, include=True):
            self.lmin, self.lmax, self.include = lmin, lmax, include

    class _Set:
        def __init__(self, masks):
            self.masks = masks

    mset = _Set([_Mask(4999.6, 5000.4)])
    res = _measurement(maskset_key=ms.mask_key(mset))
    assert not res.stale_against(mset)
    mset.masks[0].lmax = 5000.5              # an edge was dragged
    assert res.stale_against(mset)
    # a measurement made before the fingerprint existed never claims staleness
    assert not _measurement(maskset_key="").stale_against(mset)


def test_is_ion_stage():
    assert sp.is_ion_stage("Fe I", 1)
    assert sp.is_ion_stage("Fe II", 2)
    assert not sp.is_ion_stage("Fe II", 1)
    assert not sp.is_ion_stage("Fe", 1)
    assert not sp.is_ion_stage("Fe III", 2)
