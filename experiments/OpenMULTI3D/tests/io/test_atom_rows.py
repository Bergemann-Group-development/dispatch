"""The level and transition objects of a model atom (``atom.levels``, ``atom.transitions``).

A transition's ``lower`` and ``upper`` are the levels its ``i`` and ``j`` columns name, so
that ``run.plot_depart(line.lower)`` draws the lower level of the line. The atom is a
three-level table written here.
"""
import pandas as pd

from tumult.atom import Atom, Level


def three_level_atom():
    lvl = pd.DataFrame(dict(ev=[0.0, 2.1, 3.2], g=[2., 4., 2.], ion=[1, 1, 2],
                            label=['Level 1 = 3s2S', 'Level 2 = 3p2P*', 'Level 3 = 1S']))
    bbt = pd.DataFrame(dict(i=[0], j=[1], lam0=[5895.9], evi=[0.0], loggf=[-0.2]))
    return Atom.from_pandas('NA', 6.3, 23.0, lvl, bbt=bbt)


def test_a_transition_ends_on_the_levels_its_table_row_names():
    atom = three_level_atom()
    tran = atom.transitions[0]
    assert isinstance(tran.lower, Level) and isinstance(tran.upper, Level)
    assert (tran.lower.i, tran.upper.i) == (0, 1)
    assert tran.lam0 == 5895.9 and tran.loggf == -0.2
    assert tran.upper.label == 'Level 2 = 3p2P*' and tran.upper.g == 4.0
    assert tran.lower.atom is atom


def test_levels_index_like_a_list():
    levels = three_level_atom().levels
    assert len(levels) == 3 and levels[-1].ion == 2
    assert [lvl.i for lvl in levels[1:]] == [1, 2]
