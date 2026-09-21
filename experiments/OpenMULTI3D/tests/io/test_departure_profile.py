"""Reading departure coefficients back out of a solved run
(:meth:`tumult.departures.Departures.profile`).

What the reader has to get right, and what the tests below stand on:

* the saved populations are ``(nx, ny, nz, nk)`` with the LTE populations and the
  departure coefficients as the last axis of the patch record (``[..., 0]`` and
  ``[..., 1]``, the layout :meth:`tumult.run.Run.get_pop` reads), and the depth axis is
  stored as the run meshed it, which for an atmosphere is deepest first;
* a profile is against ``log tau_500`` in ascending order, so a plot of it runs from the
  surface on the left to the bottom of the atmosphere on the right;
* several columns are one profile only through their populations: the departure
  coefficient of a mean spectrum is the averaged NLTE population over the averaged LTE
  one, not the average of the per-column ratio, which a single hot column would dominate.
  The two differ whenever the LTE populations differ between columns, which is the case
  the cube tests below are built on.

The populations are stubbed: a real cube is a run output directory, which the repository
does not carry, and nothing about the arithmetic above needs one.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from tumult.departures import Departures


def _run(ltau, lte, dep):
    """A stand-in for the solved run: its saved patch record and its depth scale."""
    record = np.stack([np.asarray(lte, float), np.asarray(dep, float)], axis=-1)
    return SimpleNamespace(read_patch_save=lambda kind, zdim=None: record,
                           atmos=SimpleNamespace(ltau=np.asarray(ltau, float)))


@pytest.fixture
def solved(tmp_path):
    return Departures(folder=tmp_path, name=tmp_path.name, atmosphere="model",
                      species=("Mg",))


def _read(monkeypatch, solved, run):
    monkeypatch.setattr("tumult.run.Run.read", staticmethod(lambda folder, **kw: run))
    return solved.profile()


def test_one_column_is_read_level_by_level_against_ascending_optical_depth(
        monkeypatch, solved):
    # one column, four depth points stored deepest first, two levels
    ltau = [[[2.0, 1.0, 0.0, -1.0]]]
    lte = [[[[10.0, 20.0], [10.0, 20.0], [10.0, 20.0], [10.0, 20.0]]]]
    dep = [[[[1.0, 1.0], [1.0, 0.9], [1.0, 0.5], [1.0, 0.2]]]]

    profile = _read(monkeypatch, solved, _run(ltau, lte, dep))

    assert profile.species == "Mg" and profile.columns == 1
    assert profile.nk == 2 and profile.b.shape == (2, 4)
    # ascending optical depth: the outermost point first
    assert list(profile.ltau) == [-1.0, 0.0, 1.0, 2.0]
    assert list(profile.b[0]) == [1.0, 1.0, 1.0, 1.0]          # a level in LTE throughout
    assert list(profile.b[1]) == [0.2, 0.5, 0.9, 1.0]          # depleted towards the surface


def test_columns_are_averaged_through_their_populations(monkeypatch, solved):
    """Two columns on one depth scale: b is the ratio of the sums, not the mean ratio.

    The columns hold LTE populations 1 and 9 with departures 1 and 2, so the mean spectrum
    sees (1*1 + 9*2)/(1 + 9) = 1.9, while the mean of the two coefficients would be 1.5.
    """
    ltau = [[[1.0, 0.0], [1.0, 0.0]]]
    lte = [[[[1.0], [1.0]], [[9.0], [9.0]]]]
    dep = [[[[1.0], [1.0]], [[2.0], [2.0]]]]

    profile = _read(monkeypatch, solved, _run(ltau, lte, dep))

    assert profile.columns == 2
    assert profile.b.shape == (1, 2)
    assert np.allclose(profile.b[0], 1.9)


def test_columns_on_different_depth_scales_meet_on_the_mean_one(monkeypatch, solved):
    """Each column is interpolated onto the mean log tau scale before it is added.

    The two columns are offset by one in log tau and carry the same populations, so the
    mean scale lies halfway between them and the departure coefficient interpolated there
    is the average of each column's neighbouring points.
    """
    ltau = [[[2.0, 1.0, 0.0], [3.0, 2.0, 1.0]]]
    lte = [[[[1.0], [1.0], [1.0]], [[1.0], [1.0], [1.0]]]]
    dep = [[[[1.0], [0.6], [0.2]], [[1.0], [0.6], [0.2]]]]

    profile = _read(monkeypatch, solved, _run(ltau, lte, dep))

    assert list(profile.ltau) == [0.5, 1.5, 2.5]
    # the first column spans 0..2 and gives (0.4, 0.8, 1.0 held at its deepest point), the
    # second spans 1..3 and gives (0.2 held at its outermost point, 0.4, 0.8)
    assert np.allclose(profile.b[0], [0.3, 0.6, 0.9])


def test_a_solve_without_saved_populations_says_so(tmp_path):
    with pytest.raises(ValueError, match="no saved populations"):
        Departures(folder=tmp_path, name=tmp_path.name, atmosphere="model").profile()


def test_the_convergence_report_is_read_per_level(tmp_path):
    (tmp_path / "convergence_report_000.txt").write_text(
        "   lvl     maxconv\n     1 -8.5001E-05\n     2  1.2000E-03\n")
    (tmp_path / "convergence_report_001.txt").write_text(
        "   lvl     maxconv\n     1 -2.0000E-04\n     2  4.0000E-04\n")

    levels = Departures(folder=tmp_path, name=tmp_path.name,
                        atmosphere="model").convergence_levels()

    # the largest correction over the ranks, as a magnitude, level by level
    assert np.allclose(levels, [2.0e-4, 1.2e-3])
