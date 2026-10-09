"""Conversion of MULTI-format model atoms to the table format: atom_convert.x.

dispatch.x reads every model atom through read_pandas and converts a MULTI atom first
(HOP/atom_convert_mod.f90); ``Atom.from_old`` runs the same conversion through
``atom_convert.x``. A wrong conversion changes the atom a run uses without any error, so
these tests check each MULTI construct that has no one-to-one table column against the
value the MULTI conventions define for it: Ritz wavelengths from the level energies in
cm^-1, component positions lam0*exp(-dq*lam0), hydrogenic cross-sections alpha0 (lam/lam0)^3
on nodes linear in wavelength, and gencol's interpolation of collision tables (monotone
Hermite, constant outside the table). A wrong unit, a dropped or misplaced component, a
reversed bf table or a union grid that moves a record's own values would each fail one
test. The atom is a four-level Na-like model written here, so the tests need only the
standalone build (``make -C HOP/standalone``). The converter writes each column with the
precision ``tumult.atom.Atom.write`` gives it, so values agree to half a unit in the last
written digit: wavelengths to 5e-4 AA (through ``tumult.physics.air2vac`` for lines),
loggf to 5e-4 dex, damping and cross-sections to 5e-4 relative (three decimals in e
format) and collision data to 5e-5 relative (four decimals).
"""
import subprocess

import numpy as np
import pytest

from tumult.atom import CONVERTER, Atom
from tumult.physics import air2vac

pytestmark = pytest.mark.skipif(not CONVERTER.is_file(),
                                reason='needs HOP/standalone/build/atom_convert.x')

# levels [cm^-1]; line 2 has two components (dq [cm^-1], weight, f, damping); continuum 1
# is tabulated from its edge, continuum 2 hydrogenic; CE and CH sit on different TEMP grids
ATOM = """\
* provenance note
NA
* second note
  6.30  22.99
  4  2  2  0
      0.000   2.0  'g1'  1
  16967.630   2.0  'p1'  1
  16973.380   4.0  'p3'  1
  41449.440   1.0  'ii'  2
  2  1  0.320  21  4.0  0.0  0  6.0e+07  1.5  0.0  VOIGT
  3  1  0.640  21  4.0  0.0  2  6.2e+07  1.5  0.0  VOIGT
   -0.050  0.400  0.256  6.0e+07  1.5  0.0
    0.030  0.600  0.384  6.5e+07  1.5  0.0
  4  1  1.0e-19  3  -1.0
   2412.0  1.0e-19
   2000.0  8.0e-20
   1500.0  5.0e-20
  4  2  2.0e-18  4  1000.0
GENCOL
TEMP
  3  3000.0  5000.0  8000.0
CE
  1  2  1.0e-08  2.0e-08  3.0e-08
TEMP
  2  4000.0  6000.0
CH
  1  3  1.0e-10  2.0e-10
END
"""
ECM = np.array([0.0, 16967.630, 16973.380, 41449.440])
# air wavelengths are written to 1e-3 AA; air2vac scales the rounding by n_air
LAM_TOL = 5e-4 * 1.0003


@pytest.fixture(scope='module')
def converted(tmp_path_factory):
    path = tmp_path_factory.mktemp('atom') / 'atom.na4'
    path.write_text(ATOM)
    return Atom.from_old(path)


def test_multi_comments_are_not_copied(tmp_path):
    """The MULTI comments describe MULTI's layout and units, so the converted file
    carries only its own provenance line."""
    src, out = tmp_path / 'atom.na4', tmp_path / 'na4.txt'
    src.write_text(ATOM)
    subprocess.run([str(CONVERTER), str(src), str(out)], check=True)
    comments = [t for t in out.read_text().splitlines() if t.lstrip().startswith('*')]
    assert len(comments) == 1 and comments[0].startswith('* converted from')


def test_line_wavelength_follows_the_level_energies(converted):
    """A MULTI line sits at its Ritz wavelength: the table holds the air wavelength
    that read_pandas converts back, so its vacuum value is 1e8/(E_j - E_i)."""
    assert air2vac(converted.bbt.lam0[0]) == pytest.approx(1e8 / (ECM[1] - ECM[0]), abs=LAM_TOL)
    assert converted.bbt.loggf[0] == pytest.approx(np.log10(2.0 * 0.32), abs=5e-4)


def test_components_become_lines_behind_the_classic_ones(converted):
    """Line 2 keeps its place with its first component and the second is appended:
    wavelength lam0*exp(-dq*lam0), gf times the weight, the component's damping."""
    bbt = converted.bbt
    assert len(bbt) == 3 and list(bbt.i) == [0, 0, 0] and list(bbt.j) == [1, 2, 2]
    lam0 = 1e8 / (ECM[2] - ECM[0])
    for row, dq, wt, grad in ((1, -0.05, 0.4, 6.0e7), (2, 0.03, 0.6, 6.5e7)):
        lam_vac = lam0 * np.exp(-dq * lam0 * 1e-8)
        assert air2vac(bbt.lam0[row]) == pytest.approx(lam_vac, abs=LAM_TOL)
        assert bbt.loggf[row] == pytest.approx(np.log10(2.0 * 0.64 * wt), abs=5e-4)
        assert bbt.grad[row] == pytest.approx(grad, rel=5e-4)


def test_hydrogenic_continuum_is_tabulated_on_the_multi_nodes(converted):
    """nq nodes linear in wavelength from the edge to qmax, alpha0 (lam/lam0)^3."""
    bft = converted.bft[converted.bft.i == 1]
    edge = 1e8 / (ECM[3] - ECM[1])
    lam = np.linspace(edge, 1000.0, 4)
    np.testing.assert_allclose(bft.lam_air, lam, rtol=0, atol=5e-4)
    np.testing.assert_allclose(bft.alpha, 2.0e-18 * (lam / edge)**3, rtol=5e-4)
    tab = converted.bft[converted.bft.i == 0]
    np.testing.assert_allclose(tab.lam_air, [2412.0, 2000.0, 1500.0], rtol=0, atol=5e-4)


def test_collision_grids_merge_exactly_on_each_record_s_nodes(converted):
    """Records on different TEMP grids share the union grid. Each keeps its own values
    on its own nodes and gencol's constant continuation outside them; CH, interpolated
    linearly by MULTI convention, takes the linear value between its two nodes."""
    col = converted.col
    temps = np.array([float(t) for t in col.columns[3:]])
    np.testing.assert_array_equal(temps, [3000.0, 4000.0, 5000.0, 6000.0, 8000.0])
    ce = col[col.key == 'CE'].iloc[0, 3:].to_numpy(float)
    ch = col[col.key == 'CH'].iloc[0, 3:].to_numpy(float)
    np.testing.assert_allclose(ce[[0, 2, 4]], [1e-8, 2e-8, 3e-8], rtol=5e-5)
    np.testing.assert_allclose(ch, [1e-10, 1e-10, 1.5e-10, 2e-10, 2e-10], rtol=5e-5)
