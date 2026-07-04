"""Pin model-atom reading (AtomClass and FastAtom on data/demo/atoms/atom.O).

The partition function goes through EK and saha_const, so this golden is the canary for
the Phase 1 constants consolidation (atom.py/trans.py/level.py -> tumult.const).
"""
import numpy as np
import pytest

from cases import DEMO_ATOM, atom_cases
from conftest import assert_matches_golden, load_golden

pytestmark = pytest.mark.skipif(not DEMO_ATOM.exists(), reason="demo atom file not present")


def test_matches_golden():
    assert_matches_golden(atom_cases(), load_golden("atom"))


def test_slow_and_fast_reader_agree():
    """AtomClass and FastAtom parse the same file - level energies and line
    wavelengths must agree (both go vac->air through the same conversion)."""
    from tumult.atom import AtomClass
    from tumult.fast_atom import FastAtom

    slow = AtomClass(str(DEMO_ATOM), read_bf=True)
    fast = FastAtom(str(DEMO_ATOM))
    np.testing.assert_allclose(
        np.array([lvl.ev for lvl in slow.level]),
        np.asarray(fast.lvl.ev.values, float), rtol=1e-10)
    np.testing.assert_allclose(
        np.array([t.lam0 for t in slow.trans[: slow.nline]]),
        np.asarray(fast.bbt.lam0.values, float), rtol=1e-10)
