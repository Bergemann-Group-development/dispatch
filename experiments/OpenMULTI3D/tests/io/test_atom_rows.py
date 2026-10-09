"""The level and transition objects of a model atom (``atom.levels``, ``atom.bbtrans``,
``atom.bftrans``, ``atom.colls``).

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
    bft = pd.DataFrame(dict(i=[0, 0, 1, 1, 1], lam_air=[2412., 2000., 4083., 3500., 3000.],
                            alpha=[1e-19, 8e-20, 1e-17, 7e-18, 5e-18]))
    col = pd.DataFrame({'key': ['CE', 'CE', 'CI'], 'i': [0, 0, 1], 'j': [1, 1, 2],
                        '5000.0': [1e-9, 2e-9, 1e-12], '8000.0': [9e-10, 1e-9, 3e-12]})
    return Atom.from_pandas('NA', 6.3, 23.0, lvl, bbt=bbt, bft=bft, col=col)


def test_a_transition_ends_on_the_levels_its_table_row_names():
    atom = three_level_atom()
    tran = atom.bbtrans[0]
    assert isinstance(tran.lower, Level) and isinstance(tran.upper, Level)
    assert (tran.lower.i, tran.upper.i) == (0, 1)
    assert tran.lam0 == 5895.9 and tran.loggf == -0.2
    assert tran.upper.label == 'Level 2 = 3p2P*' and tran.upper.g == 4.0
    assert tran.lower.atom is atom


def test_levels_index_like_a_list():
    levels = three_level_atom().levels
    assert len(levels) == 3 and levels[-1].ion == 2
    assert [lvl.i for lvl in levels[1:]] == [1, 2]


def test_a_bound_free_transition_spans_the_nodes_of_one_continuum():
    atom = three_level_atom()
    bf = atom.bftrans
    assert len(bf) == atom.ncont == 2
    assert bf[1].i == 1 and list(bf[1].lam_air) == [4083., 3500., 3000.]
    assert bf[0].upper.i == bf[1].upper.i == 2   # without j: the next stage's ground level


def test_collisions_are_one_per_record():
    colls = three_level_atom().colls
    assert len(colls) == 3 and colls[1].key == 'CE'
    assert (colls[2].lower.i, colls[2].upper.i) == (1, 2)


def test_a_view_displays_and_queries_as_its_table():
    atom = three_level_atom()
    assert repr(atom.colls) == repr(atom.col)
    assert len(atom.levels.query('ion == 1')) == 2
