"""Pin the abundance numerics: curve-of-growth inversion, inverse-variance aggregation,
and the GUI-side physics that Phase 2 moves into tumult.abundance (scale conversions,
solar-table parsing, EW->abundance error propagation)."""
import numpy as np
import pytest

from cases import AGG_LINES, abundance_cases
from conftest import assert_matches_golden, load_golden


def test_matches_golden():
    assert_matches_golden(abundance_cases(), load_golden("abundance"))


def test_cog_invert_recovers_exact_grid_points():
    """The PCHIP interpolant must pass through the CoG nodes exactly."""
    from tumult.abundance import cog_invert

    a_grid = np.linspace(-0.6, 0.6, 7)
    ew_grid = 55.0 * 10 ** (0.62 * a_grid)
    for a, w in zip(a_grid, ew_grid):
        inv = cog_invert(a_grid, ew_grid, ew_meas=w)
        assert inv["a"] == pytest.approx(a, abs=1e-12)
    assert cog_invert(a_grid, ew_grid, ew_meas=ew_grid.max() * 2)["out_of_range"]
    with pytest.raises(ValueError):
        cog_invert([0.0], [10.0], ew_meas=10.0)


def test_aggregate_weighting_and_budget():
    """Inverse-variance weighting; photon+continuum shares must recompose sigma_mean."""
    from tumult.abundance import aggregate

    res = aggregate(AGG_LINES)
    inc = [ln for ln in AGG_LINES if ln.include]
    assert res["n"] == len(inc)
    w = np.array([1 / ln.sigma_a**2 for ln in inc])
    a = np.array([ln.a for ln in inc])
    assert res["a_mean"] == pytest.approx((w * a).sum() / w.sum(), rel=1e-12)
    assert res["sigma_mean"] == pytest.approx(np.sqrt(1 / w.sum()), rel=1e-12)
    # the split recomposes exactly because sigma_a**2 = phot**2 + cont**2 for these lines
    assert np.hypot(res["sigma_phot"], res["sigma_cont"]) == pytest.approx(
        res["sigma_mean"], rel=1e-6)
    empty = aggregate([])
    assert empty["n"] == 0 and np.isnan(empty["a_mean"])


def test_abundance_scale_conversions_are_exact():
    """A(X) = A_solar(X) + [Fe/H] + [X/Fe] and its inverse (GUI synth.py; moves in Phase 2).

    Pure float arithmetic - pinned exactly, and the round trip must be the identity.
    """
    import synth  # TumultStudio2/synth.py

    solar = {"Mg": 7.55, "Fe": 7.50, "O": 8.69}
    assert synth.xfe_to_abund(0.3, "Mg", -0.5, solar) == 7.35
    assert synth.abund_to_xfe(7.35, "Mg", -0.5, solar) == pytest.approx(0.3, abs=1e-14)
    # unknown element passes through unchanged
    assert synth.xfe_to_abund(0.3, "Xx", -0.5, solar) == 0.3
    # round trip through the override-string converter
    text = "Mg=0.30 O=-0.15 Xx=0.1"
    absolute = synth.convert_abundances(text, to_relative=False, feh=-0.5, solar=solar)
    back = synth.convert_abundances(absolute, to_relative=True, feh=-0.5, solar=solar)
    assert back.split() == ["Mg=0.3", "O=-0.15", "Xx=0.1"]


def test_solar_table_parser(tmp_path):
    """abund_table's parsing core (moves to tumult.abundance.read_solar_abundances)."""
    import synth

    f = tmp_path / "abund.txt"
    f.write_text("# comment\nH 12.00\nFe 7.50\nMG 7.55\nbad line\nNa not_a_number\n")
    table = synth.abund_table({"abund_file": str(f)})
    assert table == {"H": 12.00, "Fe": 7.50, "Mg": 7.55}
    assert synth.abund_table({}) == {}
    assert "Fe" not in synth.abund_elements({"abund_file": str(f)})


def test_sigma_a_from_rew_guards():
    import abundance as gui_ab

    class FakeProf:
        mAA = 50.0
        rew_err_phot = 0.02
        rew_err_cont = 0.01

    sa, sap, sac = gui_ab.sigma_a_from_rew(FakeProf(), drewda=0.4)
    assert sap == pytest.approx(0.05) and sac == pytest.approx(0.025)
    assert sa == pytest.approx(np.hypot(0.05, 0.025))
    assert gui_ab.sigma_a_from_rew(FakeProf(), drewda=0.0) == (np.inf, np.inf, np.inf)


def test_velocity_hwhm_uses_c():
    import abundance as gui_ab

    # HWHM_v = 0.5 * FWHM / lam0 * c
    assert gui_ab.velocity_hwhm(0.2, 5000.0) == pytest.approx(
        0.5 * 0.2 / 5000.0 * 299792.458)
    assert np.isnan(gui_ab.velocity_hwhm(float("nan"), 5000.0))
    assert np.isnan(gui_ab.velocity_hwhm(0.2, -1.0))


def test_screen_flags():
    import abundance as gui_ab

    flags, ok = gui_ab.screen(50.0, 0.05, 5000.0)
    assert ok and flags == ()
    assert "weak" in gui_ab.screen(1.0, 0.0, 5000.0)[0]
    assert "blended" in gui_ab.screen(50.0, 0.9, 5000.0)[0]
    assert "strong" in gui_ab.screen(5000.0, 0.0, 5000.0)[0]
