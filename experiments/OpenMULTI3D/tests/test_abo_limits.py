"""Series limits of tumult.abo from level designations, against the limits Barklem, Piskunov
& O'Mara (2000, CDS J/A+AS/142/467) used for Fe I (their table1, IP 63700 cm^-1)."""
import pandas as pd
import pytest

from tumult.abo import orbitals, parent, parent_terms, series_limit

# the Fe II terms involved, lowest level of each (NIST, cm^-1)
FE2 = pd.DataFrame([('3d6.(5D).4s', 'a 6D', '9/2', 0.0), ('3d7', 'a 4F', '9/2', 1872.5998),
                    ('3d6.(5D).4s', 'a 4D', '7/2', 7955.3186), ('3d7', 'a 4P', '5/2', 13474.4474),
                    ('3d6.(3P2).4s', 'b 4P', '5/2', 20830.5534), ('3d6.(3P1).4s', 'c 4P', '5/2', 49100.956),
                    ('3d5.4s2', 'a 6S', '5/2', 23317.6351)], columns=['config', 'term', 'J', 'ecm'])
MG2 = pd.DataFrame([('2p6.3s', '2S', '1/2', 0.0), ('2p6.3p', '2P*', '1/2', 35669.31)],
                   columns=['config', 'term', 'J', 'ecm'])


@pytest.mark.parametrize('config, limit', [
    ('3d6.4s2', 63700.0),                        # a5D: 4s, open parent term -> a6D
    ('3d6.(5D).4s.4p.(3P*)', 63700.0),           # z levels: 4p on 3d6(5D)4s -> a6D
    ('3d7.(4F).4s', 65572.0),                    # a5F: a4F
    ('3d7.(4F).4p', 65572.0),
    ('3d6.(5D).4s\\ (4D).5s', 71655.0),          # VALD's '\ (4D)': a4D
    ('3d7.(4P).4s', 77174.0),
    ('3d5.(6S).4s2.4p', 87017.0),                # a6S, NIST without the core term
    ('3d6.(a3P).4s.4p.(3P*)', 84530.0),          # 'a3P' is NIST's lower 3P2 (b4P); BPM took b2P
])
def test_fe1_limits_follow_bpm(config, limit):
    lim, _ = series_limit(parent(orbitals(config)), parent_terms(FE2), 63700.0)
    assert lim == pytest.approx(limit, abs=1.0)


def test_closed_shells_are_ignored():
    """NIST writes Mg II's 2p6, Mg I designations leave it out."""
    lim, n = series_limit(parent(orbitals('3s.4s')), parent_terms(MG2), 61671.05)
    assert lim == pytest.approx(61671.05) and n == 1
