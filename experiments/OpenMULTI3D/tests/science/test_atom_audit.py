"""compare_linelist attributes a line-list line to atom levels only through levels of a
compatible LS term.

Mg I 3p2 3P (J = 0, 1, 2 at 7.168-7.176 eV) lies within 3 meV of 3s7s 3S (7.173 eV). An
atom without 3p2 3P must not receive the 3p 3P* - 3p2 3P multiplet (log gf up to +0.4)
as part of its 3p 3P* - 7s 3S transition (log gf about -2): the comparison then reports
a 2.4 dex gf deficit that does not exist. Energies are NIST ASD values; the line-list
rows are VALD's 2779.83 and 2778.14 AA lines. The synthetic atom keeps only the levels
the case needs.
"""
import numpy as np
import pandas as pd

from tumult.atom_audit import HC_AA_EV, compare_linelist
from tumult.fast_atom import PandasAtom
from tumult.physics import vac2air


def mg_atom_without_3p2():
    lvl = pd.DataFrame(dict(ev=[0.0, 2.713916, 7.173192, 7.646299], g=[1., 9., 3., 2.], ion=[1, 1, 1, 2],
                            label=['3s1S', '3p3P*', '7s3S', '3s2S']))
    lam = float(vac2air(HC_AA_EV / (lvl.ev[2] - lvl.ev[1])))
    bbt = pd.DataFrame(dict(i=[1], j=[2], lam0=[lam], evi=[lvl.ev[1]], loggf=[-1.6]))
    return PandasAtom.from_pandas('MG', 7.5, 24.3, lvl, bbt=bbt)


def nist_levels():
    rows = [('3s.3p', '3P*', 0, 2.709127), ('3s.3p', '3P*', 1, 2.711614), ('3s.3p', '3P*', 2, 2.716662),
            ('3p2', '3P', 0, 7.167929), ('3p2', '3P', 1, 7.170487), ('3s.7s', '3S', 1, 7.173192),
            ('3p2', '3P', 2, 7.175513)]
    tab = pd.DataFrame(rows, columns=['config', 'term', 'J', 'ev'])
    tab['g'] = 2 * tab.J + 1
    tab['ion'] = 1
    return tab


def test_line_to_a_term_the_atom_lacks_is_not_given_to_a_neighbour():
    llist = pd.DataFrame(dict(lam0=[2779.831, 2778.135], evi=[2.7166, 2.7116], loggf=[0.375, -2.23]))
    comp, _ = compare_linelist(mg_atom_without_3p2(), llist, 1, nist_levels=nist_levels())
    row = comp.loc[(1, 2)]
    assert row.n_components == 1
    np.testing.assert_allclose(np.log10(row.gf_list), -2.23, atol=1e-9)


def test_line_terms_choose_between_degenerate_singlet_and_triplet():
    """Mg I 4f 1F* and 4f 3F* lie 0.03 meV apart, far below the energy a wavelength and a
    lower level give, so energy alone hands 3d 3D - 4f 3F* components to whichever NIST
    level is nearest. With the line list's terms the component lands on the atom's 4f 3F*.
    """
    lvl = pd.DataFrame(dict(ev=[0.0, 5.945964, 6.779070, 6.779104, 7.646299], g=[1., 15., 7., 21., 2.],
                            ion=[1, 1, 1, 1, 2], label=['3s1S', '3d3D', '4f1F*', '4f3F*', '3s2S']))
    lam = float(vac2air(HC_AA_EV / (lvl.ev[3] - lvl.ev[1])))
    bbt = pd.DataFrame(dict(i=[1], j=[3], lam0=[lam], evi=[lvl.ev[1]], loggf=[1.07]))
    atom = PandasAtom.from_pandas('MG', 7.5, 24.3, lvl, bbt=bbt)
    nist = pd.DataFrame(dict(config='3s.4f', term=['1F*', '3F*', '3F*', '3F*'], J=[3, 2, 3, 4],
                             ev=[6.779070, 6.779104, 6.779104, 6.779104]))
    nist['g'] = 2 * nist.J + 1
    nist['ion'] = 1
    # a component whose upper energy, from its rounded wavelength, falls on the 1F* level
    lam_line = float(vac2air(HC_AA_EV / (6.779070 - lvl.ev[1])))
    llist = pd.DataFrame(dict(lam0=[lam_line], evi=[lvl.ev[1]], loggf=[0.7], term_lo=['3D'], term_up=['3F*']))
    comp, _ = compare_linelist(atom, llist, 1, nist_levels=nist)
    assert (1, 3) in comp.index and comp.loc[(1, 3)].n_components == 1
    assert (1, 2) not in comp.index
