"""Pin the window-measuring core shared by the abundance and ionisation-balance workflows:
what a measured window yields, how the dominant transition in it is identified (and when that
identification is refused), and that the masks themselves come back untouched.
"""
import numpy as np
import pandas as pd
import pytest

from tumult import linediag
from tumult.linemask import LineMask
from tumult import Spectrum

LAM1, LAM2 = 5000.0, 5010.0


def _spectrum(depths=(0.45, 0.30)):
    x = np.arange(4998.0, 5012.0, 0.005)
    y = np.ones_like(x)
    for lam, d in zip((LAM1, LAM2), depths):
        y -= d * np.exp(-0.5 * ((x - lam) / 0.08) ** 2)
    return Spectrum(x, y)


def _catalog(rows):
    """``rows`` of (ion, lam0, evi, loggf, strength)."""
    return pd.DataFrame(
        [dict(zip(("ion", "lam0", "evi", "loggf", "strength"), r)) for r in rows])


def _masks():
    return [LineMask(lam0=lam, lmin=lam - 0.4, lmax=lam + 0.4,
                     smin=lam - 0.9, smax=lam + 0.9) for lam in (LAM1, LAM2)]


# --- identification -------------------------------------------------------------------------

def test_identify_takes_the_strongest_transition_in_the_window():
    cat = _catalog([("Fe I", 4999.8, 1.0, -2.0, 0.1),
                    ("Fe II", 5000.1, 3.0, -1.0, 9.0),
                    ("Ni I", 5000.3, 2.0, -1.5, 0.5)])
    ident, flags = linediag._identify(cat, 4999.6, 5000.4)
    assert ident["ion"] == "Fe II" and ident["e_low"] == 3.0 and flags == ()


def test_identify_flags_a_window_another_species_dominates():
    cat = _catalog([("Fe I", 4999.9, 1.0, -2.0, 0.2), ("Ni I", 5000.1, 2.0, -1.0, 5.0)])
    ident, flags = linediag._identify(cat, 4999.6, 5000.4, element="Fe")
    assert ident["ion"] == "Ni I" and flags == ("off_element",)
    # without an element to check against, the same window is simply identified
    assert linediag._identify(cat, 4999.6, 5000.4)[1] == ()
    # and a window the element does dominate is never flagged
    assert linediag._identify(cat, 4999.6, 5000.0, element="Fe")[1] == ()


def test_identify_reports_an_empty_window():
    cat = _catalog([("Fe I", 5100.0, 1.0, -2.0, 1.0)])
    ident, flags = linediag._identify(cat, 4999.6, 5000.4, element="Fe")
    assert ident == {} and flags == ("no_lines",)
    # no catalog at all is not a failure: the EW is still measurable, the identity unknown
    assert linediag._identify(None, 4999.6, 5000.4) == ({}, ())


# --- measurement ------------------------------------------------------------------------------

def test_measure_masks_measures_each_window():
    masks = _masks()
    cat = _catalog([("Fe I", LAM1, 1.0, -1.0, 1.0), ("Fe II", LAM2, 3.0, -1.5, 1.0)])
    lines = linediag.measure_masks(_spectrum(), masks, profile="gauss", snr=200.0,
                                   cont_pct=0.5, catalog=cat, element="Fe")
    # 0.45 * sqrt(2 pi) * 0.08 AA = 90.2 maa, and 60.2 maa for the shallower line
    assert lines[0].ew == pytest.approx(90.2, abs=0.5)
    assert lines[1].ew == pytest.approx(60.2, abs=0.5)
    assert [ln.ion for ln in lines] == ["Fe I", "Fe II"]
    assert [ln.e_low for ln in lines] == [1.0, 3.0]
    for ln in lines:
        assert ln.fwhm == pytest.approx(2.3548 * 0.08, abs=1e-3)
        assert 0 < ln.ew_err < ln.ew
        assert ln.flags == ()


def test_measure_masks_flags_an_unfittable_window():
    m = LineMask(lam0=float("nan"), lmin=float("nan"), lmax=float("nan"))
    lines = linediag.measure_masks(_spectrum(), [m], profile="gauss", snr=None, cont_pct=0.0)
    assert lines[0].flags == ("no_window",)


def test_measure_masks_leaves_the_masks_alone():
    """The set is the human's windows: a measurement copies them and writes nothing back."""
    masks = _masks()
    before = [m.to_dict() for m in masks]
    lines = linediag.measure_masks(_spectrum(), masks, profile="gauss", snr=200.0,
                                   cont_pct=0.5)
    assert [m.to_dict() for m in masks] == before
    # each line carries its own copy of the window, so editing the set leaves it alone
    ln = lines[0]
    assert ln.fit_window == (masks[0].lmin, masks[0].lmax)
    assert ln.window == (masks[0].smin, masks[0].smax)
    masks[0].lmax = 5001.0
    assert ln.fit_window[1] == pytest.approx(5000.4)


# --- generate_mask_set: the headless "make me a diagnostic set for element X" call ---------

FE_LAMS = (5000.30, 5000.80, 5001.40)
TI_LAM = 5000.55


def _mini_master(tmp_path):
    """A two-species line list in tumult format: three Fe I lines and one Ti I."""
    rows = "".join(f"  {l:.4f}  1.011  -2.500  1.0  8.0  -7.0  1  2\n" for l in FE_LAMS)
    p = tmp_path / "mini_master.txt"
    p.write_text(f"'Fe I'\n'26.00'  {len(FE_LAMS)}\n{rows}"
                 f"'Ti I'\n'22.00'  1\n  {TI_LAM:.4f}  0.900  -1.900  1.0  8.0  -7.0  1  2\n")
    return p


def _blended_spectrum():
    x = np.linspace(4999.8, 5001.8, 3000)
    y = np.ones_like(x)
    for lam in FE_LAMS + (TI_LAM,):
        y -= 0.35 * np.exp(-0.5 * ((x - lam) / 0.05) ** 2)
    return Spectrum(x, y)


def test_spectrum_range_is_ascending_angstrom():
    lo, hi = linediag.spectrum_range_aa(_blended_spectrum())
    assert lo < hi
    assert (lo, hi) == pytest.approx((4999.8, 5001.8))


def test_generate_mask_set_selects_only_the_asked_element(tmp_path):
    """A Ti line sitting between the Fe ones must not enter an Fe set."""
    mset = linediag.generate_mask_set(
        _blended_spectrum(), "Fe", linelist=_mini_master(tmp_path),
        reference_name="synthetic", profile="gauss", ew_min=1.0, snr=100.0, cont_pct=0.5)
    assert mset.element == "Fe"
    assert [pytest.approx(m.lam0, abs=1e-6) for m in mset.masks] == list(FE_LAMS)
    assert mset.label.startswith("Fe · auto")


def test_generate_mask_set_rejects_a_non_element(tmp_path):
    with pytest.raises(ValueError, match="not an element symbol"):
        linediag.generate_mask_set(_blended_spectrum(), "Xx",
                                   linelist=_mini_master(tmp_path))


def test_generate_mask_set_reports_an_element_with_no_lines(tmp_path):
    with pytest.raises(ValueError, match="No Ca lines"):
        linediag.generate_mask_set(_blended_spectrum(), "Ca",
                                   linelist=_mini_master(tmp_path))


def test_generate_mask_set_reports_when_everything_is_too_weak(tmp_path):
    with pytest.raises(ValueError, match="above"):
        linediag.generate_mask_set(_blended_spectrum(), "Fe",
                                   linelist=_mini_master(tmp_path), ew_min=1e6)


def test_strip_ion_stage_leaves_molecules_alone():
    """Dropping the stage gives the element for an atom, but a molecule carries no stage —
    it is *not* a decomposition into elements (that is isotopes.species_elements)."""
    from tumult.isotopes import species_elements
    from tumult.linelist import strip_ion_stage, strip_ion_stages

    assert strip_ion_stage("Fe I") == "Fe"
    assert strip_ion_stage("Fe  II") == "Fe"      # whitespace normalised
    assert strip_ion_stage("CH") == "CH"          # a molecule, not an element
    assert strip_ion_stage("") == ""
    assert species_elements("CH") == ["C", "H"]   # the decomposition, for contrast
    assert list(strip_ion_stages(pd.Series(["Fe I", "Ti II", "CO"]))) == ["Fe", "Ti", "CO"]


def test_select_candidates_excludes_molecular_species(tmp_path):
    """A molecular species never matches an element symbol, so it cannot enter an Fe set."""
    from tumult.linelist import load_master

    rows = "".join(f"  {l:.4f}  1.011  -2.500  1.0  8.0  -7.0  1  2\n" for l in FE_LAMS)
    p = tmp_path / "with_molecule.txt"
    p.write_text(f"'Fe I'\n'26.00'  {len(FE_LAMS)}\n{rows}"
                 f"'CH'\n'106.00'  1\n  5000.55  0.900  -1.900  1.0  8.0  -7.0  1  2\n")
    master = load_master(str(p))
    assert len(linediag.select_candidates(master, "Fe", 4999.0, 5002.0)) == len(FE_LAMS)
    assert len(linediag.select_candidates(master, "C", 4999.0, 5002.0)) == 0
