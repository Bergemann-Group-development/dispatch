"""Diagnostic-aware model-atom pruning: the relevance measure and the pruned atom."""
import contextlib
import io
from pathlib import Path

import numpy as np

from tumult.atom_build import prune
from tumult.atom_rates import diagnostic_relevance, removable
from tumult.fast_atom import PandasAtom

ATOM_FILE = Path(__file__).resolve().parents[2] / "input" / "atoms" / "atom.o41f"


def cycle_rates(transitions, nk):
    """Rates of one depth with the given net radiative rates [(i, j, net)], i < j, no
    collisions and b = 2 everywhere."""
    nstar = np.ones((nk, 1))
    pop = 2 * nstar
    irad = np.array([t[0] for t in transitions])
    jrad = np.array([t[1] for t in transitions])
    net = np.array([t[2] for t in transitions], float)
    rij = np.where(net > 0, net, 0)[:, None] / pop[irad]
    rji = np.where(net < 0, -net, 0)[:, None] / pop[jrad]
    return dict(pop=pop, nstar=nstar, rij=rij, rji=rji, irad=irad, jrad=jrad,
                cup=np.zeros((nk * (nk - 1) // 2, 1)))


def test_influence_follows_the_flow_into_the_diagnostic_levels():
    # In statistical equilibrium, a strong cycle 0 -> 1 -> 2 -> 0 through the diagnostic
    # level 0, a branch carrying 1% of level 2's flow through 3 and 4, and a separate
    # cycle between 5 and 6
    rates = cycle_rates([(0, 1, 1.0), (1, 2, 1.0), (0, 2, -1.0),
                         (2, 3, 0.01), (3, 4, 0.01), (2, 4, -0.01),
                         (5, 6, 1.0), (5, 6, -1.0)], 7)
    sources = np.zeros((7, 1), bool)
    sources[0] = True
    influence, ends = diagnostic_relevance(rates, sources)
    relevance = ends.max(1)
    branch = 0.01 / 1.01
    assert np.allclose(influence[:3], 1)
    assert np.allclose(influence[3:5], branch)
    assert np.allclose(influence[5:], 0)
    assert np.allclose(relevance[:3], 1)
    # the branch transitions dominate their own levels but reach level 0 only through
    # the 1% branch
    assert np.allclose(relevance[3:6], branch)
    assert np.allclose(relevance[6:], 0)


def test_relevance_ignores_depths_without_diagnostic_levels():
    rates = cycle_rates([(0, 1, 1.0), (1, 2, 1.0), (0, 2, -1.0)], 3)
    influence, ends = diagnostic_relevance(rates, np.zeros((3, 1), bool))
    assert not influence.any() and not ends.any()


def test_pruned_atom_drops_a_level_with_its_transitions():
    with contextlib.redirect_stdout(io.StringIO()):
        atom = PandasAtom.from_old(ATOM_FILE)
    # the non-ground level with the most lines
    ends = np.concatenate([atom.bbt.i, atom.bbt.j])
    counts = np.bincount(ends, minlength=len(atom.lvl))
    counts[atom.lvl.ion.ne(atom.lvl.ion.shift()).to_numpy()] = 0
    gone = int(counts.argmax())
    keep = np.ones(len(atom.lvl), bool)
    keep[gone] = False
    out = prune(atom, np.ones(len(atom.bbt), bool), keep)

    assert len(out.lvl) == len(atom.lvl) - 1
    assert len(out.bbt) == len(atom.bbt) - counts[gone]
    assert len(out.bft) == len(atom.bft) - (atom.bft.i == gone).sum()
    assert len(out.col) == len(atom.col) - ((atom.col.i == gone) | (atom.col.j == gone)).sum()
    # the remaining transitions join the same levels as before
    old_ev = atom.lvl.ev.to_numpy()
    kept = atom.bbt[(atom.bbt.i != gone) & (atom.bbt.j != gone)]
    assert np.allclose(np.sort(out.lvl.ev.to_numpy()[out.bbt.j] - out.lvl.ev.to_numpy()[out.bbt.i]),
                       np.sort(old_ev[kept.j] - old_ev[kept.i]))
    assert out.bbt.i.max() < len(out.lvl) and out.col[['i', 'j']].to_numpy().max() < len(out.lvl)


def test_removal_budget_bounds_the_sum_per_level():
    # four transitions into level 0, each carrying 0.4 of the budget there
    ends = np.array([[0.04, 0.0], [0.04, 0.0], [0.04, 0.0], [0.0, 0.5]])
    lower = np.array([0, 0, 0, 1])
    upper = np.array([2, 3, 4, 0])
    gone = removable(ends, lower, upper, 5, budget=0.1)
    assert gone.tolist() == [True, True, False, False]
