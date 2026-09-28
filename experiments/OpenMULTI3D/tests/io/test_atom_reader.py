"""Model-atom reading by AtomClass and FastAtom, on the tracked ``input/atoms/atom.o41f``.

The reference ``tests/golden/atom.npz`` pins counts, level energies and weights, line
wavelengths, B_ij and two partition-function values. The partition function goes through
EK and saha_const, so a change to the shared constants shows up here. Reading is pure
float64 parsing and arithmetic, so the tolerance is roundoff (rtol 1e-12).
"""
import numpy as np

from tests.io.atom_reference import ATOM_FILE, atom_reference
from tests.support.golden import assert_matches_golden, load_golden


def test_matches_reference():
    assert_matches_golden(atom_reference(), load_golden("atom"), rtol=1e-12)


def test_slow_and_fast_reader_agree():
    """AtomClass and FastAtom parse the same file - level energies and line
    wavelengths must agree (both go vac->air through the same conversion)."""
    from tumult.atom import AtomClass
    from tumult.fast_atom import FastAtom

    slow = AtomClass(str(ATOM_FILE), read_bf=True)
    fast = FastAtom(str(ATOM_FILE))
    np.testing.assert_allclose(
        np.array([lvl.ev for lvl in slow.level]),
        np.asarray(fast.lvl.ev.values, float), rtol=1e-10)
    np.testing.assert_allclose(
        np.array([t.lam0 for t in slow.trans[: slow.nline]]),
        np.asarray(fast.bbt.lam0.values, float), rtol=1e-10)
