"""The parts of tumult.line_quadrature that need no HOP: the thinning invariant and the
atom copy, whose contract is that only the ``nq qmax q0`` columns change."""
from pathlib import Path

import numpy as np

from tumult.line_quadrature import (read_classic_atom, thin_linear, write_atom_copy,
                                    write_atom_linelist)
from tumult.linelist import LineList

ATOM = Path(__file__).resolve().parents[2] / "input" / "atoms" / "atom.ba06"


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


def test_atom_linelist_reads_back(tmp_path):
    atom = read_classic_atom(ATOM)
    ll = LineList.from_tumult_format(write_atom_linelist(atom, tmp_path / "ll.txt"))
    assert len(ll) == atom.lines.size
    assert set(ll.index.get_level_values("ion")) == {"Ba II"}
    # Ba II 4554: resonance line, log gf = 0.14 in the atom's comment
    row = ll.iloc[np.argmin(np.abs(ll["lam0"].values - 4554.03))]
    assert abs(row["evi"]) < 1e-6 and abs(row["loggf"] - 0.14) < 5e-3
