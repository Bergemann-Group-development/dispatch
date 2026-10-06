"""tumult.hop.read_atmosphere: dispatch.x's atmosphere reader through libhop.

The requirement is that Python sees the model exactly as the Fortran reader holds it.
For a MARCS model the expected values are the file's own columns, parsed here from the
text independently of the Fortran.
"""
from pathlib import Path

import numpy as np
import pytest

from tests.native.hop_reference import HOP_LIB

pytestmark = pytest.mark.skipif(
    not HOP_LIB.exists(),
    reason="libhop.so not built - run 'make hop_standalone'")

MARCS = (Path(__file__).resolve().parents[2] / 'input' / 'atmos' /
         'p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod')


def marcs_columns(path):
    lines = path.read_text().splitlines()
    i = next(j for j, line in enumerate(lines) if 'Number of depth points' in line)
    nz = int(lines[i].split()[0])
    first = np.array([lines[i + 3 + k].split() for k in range(nz)], float)
    second = np.array([lines[i + 4 + nz + k].split() for k in range(nz)], float)
    return first[:, 4], first[:, 5], second[:, 3]          # T, Pe, rho


def test_marcs_columns_match_the_file():
    from tumult.hop import read_atmosphere
    atm = read_atmosphere(MARCS)                            # format from the name
    temp, pe, rho = marcs_columns(MARCS)
    assert atm.shape == (1, 1, temp.size)
    assert atm.has_pe and atm.has_rho and atm.feh == 0.0 and atm.abund.size == 92
    for got, want in ((atm.ltemp, temp), (atm.lpe, pe), (atm.lrho, rho)):
        np.testing.assert_allclose(np.exp(got[0, 0].astype(float)), want, rtol=1e-6)


def test_missing_file_raises_with_the_reader_message():
    from tumult.hop import HopError, read_atmosphere
    with pytest.raises(HopError, match='File not found'):
        read_atmosphere('no_such_model.mod')
