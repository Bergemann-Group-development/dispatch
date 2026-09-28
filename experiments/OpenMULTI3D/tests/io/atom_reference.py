"""The model-atom reading reference: AtomClass and FastAtom on the tracked O I/II atom."""
from pathlib import Path

import numpy as np

ATOM_FILE = Path(__file__).resolve().parents[2] / "input" / "atoms" / "atom.o41f"


def atom_reference():
    """Counts, level energies [eV] and weights, line wavelengths [AA, air], B_ij and two
    partition-function values from both readers."""
    from tumult.atom import AtomClass
    from tumult.fast_atom import FastAtom

    atom = AtomClass(str(ATOM_FILE), read_bf=True)
    out = {
        "nk": float(atom.nk),
        "nline": float(atom.nline),
        "ncont": float(atom.ncont),
        "abnd": atom.abnd,
        "awgt": atom.awgt,
        "level_ev": np.array([lvl.ev for lvl in atom.level]),
        "level_g": np.array([lvl.g for lvl in atom.level]),
        "line_lam0": np.array([t.lam0 for t in atom.trans[: atom.nline]]),
        "line_bij": np.array([t.bij for t in atom.trans[: atom.nline]]),
        # the partition function goes through EK and saha_const
        "pf_5000_1e12": atom.partitionFunction(5000.0, 1e12),
        "pf_8000_1e13": atom.partitionFunction(8000.0, 1e13),
    }
    fast = FastAtom(str(ATOM_FILE))
    out["fast_lam0"] = np.asarray(fast.bbt.lam0.values, dtype=float)
    out["fast_ev"] = np.asarray(fast.lvl.ev.values, dtype=float)
    return out
