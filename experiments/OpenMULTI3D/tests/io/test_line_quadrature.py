"""The parts of tumult.line_quadrature that need no HOP opacities: the thinning invariant,
the atom copy, whose contract is that only the ``nq qmax q0`` columns change, the lines the
generator evaluates, and the grid several atoms share. Reading the lines needs the built
``atom_convert.x``."""
from pathlib import Path

import numpy as np
import pytest

from tumult.line_quadrature import (PreparedLineQuadrature, combine_line_quadratures,
                                    read_classic_atom, read_components, thin_linear,
                                    write_atom_copy, write_atom_linelist)
from tumult.atom import CONVERTER
from tumult.linelist import LineList

ATOM = Path(__file__).resolve().parents[2] / "input" / "atoms" / "atom.ba06"
HFS_ATOM = Path(__file__).resolve().parents[2] / "input" / "atoms" / "public" / "atom.mn281kbc"
needs_converter = pytest.mark.skipif(not CONVERTER.is_file(), reason="needs atom_convert.x")


def test_thin_linear_recovers_within_tolerance():
    x = np.linspace(-30, 30, 4001)
    g = np.log10(np.exp(-x**2) + 1e-2 / (1 + x**2) + 1e-4)   # Voigt-like, floored
    tol = 0.05
    keep = thin_linear(x, g, tol)
    assert keep[0] == 0 and keep[-1] == x.size - 1
    err = np.abs(np.interp(x, x[keep], g[keep]) - g)
    assert err.max() <= tol + 1e-12
    assert 10 < keep.size < 100          # neither everything nor nothing


def test_atom_copy_changes_only_the_quadrature_columns(tmp_path):
    atom = read_classic_atom(ATOM)
    half = np.full(atom.lines.size, 120.0)
    half[1] = np.nan                                  # an untouched line
    half[3] = 106.087                                 # qmax 10.6087 rounds up
    copy = write_atom_copy(atom, tmp_path / "atom.ba", half, qnorm=10.0)
    orig = ATOM.read_text().splitlines()
    new = copy.read_text().splitlines()
    assert len(orig) == len(new)
    for k, (a, b) in enumerate(zip(orig, new)):
        if k not in atom.records or k == list(atom.records)[1]:
            assert a == b
            continue
        wa, wb = a.split(), b.split()
        q = "10.7" if k == list(atom.records)[3] else "12.0"
        assert wb[3:6] == ["3", q, q]
        assert wa[:3] == wb[:3] and wa[6:] == wb[6:]
    again = read_classic_atom(copy)
    assert np.all(again.lines["qmax"][[0, 2]] == 12.0) and again.lines["qmax"][3] == 10.7
    assert again.lines["qmax"][1] == atom.lines["qmax"][1]


@needs_converter
def test_atom_linelist_reads_back(tmp_path):
    atom = read_classic_atom(ATOM)
    lines = read_components(ATOM, atom)
    ll = LineList.from_tumult_format(write_atom_linelist(atom, lines, tmp_path / "ll.txt"))
    assert len(ll) == atom.lines.size
    assert set(ll.index.get_level_values("ion")) == {"Ba II"}
    # Ba II 4554: resonance line, log gf = 0.14 in the atom's comment
    row = ll.iloc[np.argmin(np.abs(ll["lam0"].values - 4554.03))]
    assert abs(row["evi"]) < 1e-6 and abs(row["loggf"] - 0.14) < 5e-3


@needs_converter
def test_components_are_lines_of_their_own(tmp_path):
    atom = read_classic_atom(HFS_ATOM)
    lines = read_components(HFS_ATOM, atom)
    split = np.flatnonzero(atom.lines["ncomp"] > 1)
    assert split.size > 0 and len(lines) == atom.lines["ncomp"].sum()
    n = split[0]
    comps = lines[lines["parent"] == n]
    assert len(comps) == atom.lines["ncomp"][n] and comps["lam0"].nunique() > 1
    # the components share their line's total gf
    gf = atom.g[atom.lines["i"][n]] * atom.lines["f"][n]
    assert abs(np.log10((10.0 ** comps["loggf"]).sum()) - np.log10(gf)) < 0.01
    # the copy keeps the component records after their line
    half = np.full(atom.lines.size, 50.0)
    copy = read_classic_atom(write_atom_copy(atom, tmp_path / "atom.mn", half, qnorm=10.0))
    assert np.array_equal(copy.lines["ncomp"], atom.lines["ncomp"])
    assert np.all(copy.lines["qmax"] == 5.0)


def test_several_atoms_share_the_union_of_their_grids(tmp_path):
    grids = {"ba": [4554.0, 4554.1, 5853.6], "fe": [4554.05, 5000.0]}
    prepared = []
    for name, grid in grids.items():
        np.savetxt(tmp_path / f"{name}.txt", grid, header="grid")
        prepared.append(PreparedLineQuadrature(str(tmp_path), str(tmp_path / f"{name}.txt"),
                                               f"atom.{name}", f"src.{name}", "sun.mod", {}))
    combined = combine_line_quadratures(tmp_path, prepared)
    assert combined.atom_files == ["atom.ba", "atom.fe"]
    union = np.loadtxt(combined.grid_file, comments="#")
    np.testing.assert_array_equal(union, np.unique(np.concatenate(list(grids.values()))))
    assert combine_line_quadratures(tmp_path, prepared[:1]).grid_file == prepared[0].grid_file
