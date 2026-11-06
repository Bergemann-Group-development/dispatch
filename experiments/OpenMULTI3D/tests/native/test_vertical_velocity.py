"""A vertical velocity in a 1D text atmosphere shifts the line profiles of every ray.

The optional sixth column of the 'text' format is vz [km/s], positive outward. A ray of
direction cosine mu then sees the profiles shifted by mu vz, which makes the 1D line
opacity depend on the direction. The cases are the solar MARCS model written as a text
table (log tau_500 depth column, vmic 1 km/s), with a constant vz, where the answers
follow from Galilean invariance and need no reference run:

- the reader stores the column in cm/s and leaves vz zero without it;
- an all-zero column takes the static path, and the spectrum is identical;
- the emergent intensity of each ray is the static one shifted by mu vz, a blueshift for
  vz > 0, which catches a wrong sign or a shift applied to the wrong hemisphere;
- the NLTE departures do not change.

Needs the built ``dispatch.x`` and ``libhop``, the solar MARCS model, ``SIU_masterline_Bergemann2026.txt``
and ``atom.mg86g``.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

import tumult
from tests.native.hop_lib import HOP_LIB

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = ROOT / 'input' / 'atmos' / ('p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_'
                                     'o+0.00_r+0.00_s+0.00.mod')
LINELIST = Path('input/linelists/SIU_masterline_Bergemann2026.txt')
ATOM = Path('input/atoms/atom.mg86g')
C_KMS = 299792.458

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and HOP_LIB.is_file() and MARCS.is_file()
             and all((ROOT / p).is_file() for p in (LINELIST, ATOM))),
        reason='needs dispatch.x, libhop, the solar MARCS model, SIU_masterline_Bergemann2026.txt and atom.mg86g'),
]


def _write_text_model(path, vz=None):
    """The MARCS model's log tau_500, T, Pe and rho as a 'text' table, with vz [km/s]."""
    lines = MARCS.read_text().splitlines()
    i = next(j for j, line in enumerate(lines) if 'Number of depth points' in line)
    nz = int(lines[i].split()[0])
    first = np.array([lines[i + 3 + k].split() for k in range(nz)], float)
    second = np.array([lines[i + 4 + nz + k].split() for k in range(nz)], float)
    rows = ['sun_marcs', f'{nz} 0.0', '* lgtau5 T Pe rho vmic' + ('' if vz is None else ' vz')]
    for k in range(nz):
        row = f'{first[k, 2]:.4f} {first[k, 4]:.2f} {first[k, 5]:.5e} {second[k, 3]:.5e} 1.0'
        rows.append(row + ('' if vz is None else f' {vz:.4f}'))
    path.write_text('\n'.join(rows) + '\n')
    return path


def _run(tmp_path, name, vz, groups):
    work = tmp_path / name
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    _write_text_model(work / 'model.atmos', vz)
    (work / 'r.nml').write_text(f"""
&io_params           datadir='data' do_trace=F /
&atmos_params        atmos_format='text' atmos_file='model.atmos' FeH=0.0 use_ne=.false. nz=64 /
&run_params          n_windows=4 /
{groups}""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='4'), timeout=300)
    assert result.returncode == 0 and 'OpenMULTI3D terminated successfully' in result.stdout, \
        result.stdout[-2000:]
    return tumult.read(work / 'data' / 'r', verbose=0)


# Mg I b2 with its blends, sampled at 5 mA (0.29 km/s); five Gauss rays
LTE = f"""&linelist_params     linelist_file='{LINELIST}' /
&spectrum_params     daa=0.005 aa_blue=5171.0 aa_red=5175.0 /
&rt_iteration_params maxiter=0 /
&rt_output_params    angle_quad='gauss' ntheta=5 /
"""


def test_reader_stores_the_vz_column_in_cm_per_s(tmp_path):
    from hopacity import read_atmosphere
    moving = read_atmosphere(_write_text_model(tmp_path / 'v.atmos', 1.5), atmos_format='text')
    static = read_atmosphere(_write_text_model(tmp_path / 's.atmos'), atmos_format='text')
    np.testing.assert_allclose(moving.vz, 1.5e5, rtol=1e-6)
    assert np.all(static.vz == 0)
    for key in ('ltemp', 'lpe', 'lrho', 'vmic', 'zz'):
        np.testing.assert_array_equal(getattr(moving, key), getattr(static, key))


def test_zero_vz_column_gives_the_static_spectrum(tmp_path):
    """Discrete identity: the all-zero column must select the static code path."""
    static = _run(tmp_path, 'static', None, LTE)
    zero = _run(tmp_path, 'zero', 0.0, LTE)
    np.testing.assert_array_equal(np.asarray(zero.ie), np.asarray(static.ie))


def test_constant_vz_shifts_each_ray_by_mu_vz(tmp_path):
    """I_mu(lam) = I_static,mu(lam (1 + mu vz / c)).

    The static intensities are interpolated by a cubic spline on their 5 mA grid. The
    residual of that interpolation is about 2e-4 of the continuum intensity, while the
    shift of 2 km/s moves the intensity by up to 0.5 of it at mu = 1 and an opposite sign
    by more, so 1e-3 separates the two.
    """
    vz = 2.0
    static = _run(tmp_path, 'static', None, LTE)
    moving = _run(tmp_path, 'moving', vz, LTE)
    lam = np.asarray(static.lam_vac, float)
    np.testing.assert_array_equal(np.asarray(moving.lam_vac, float), lam)
    i_static, i_moving = np.asarray(static.ie, float), np.asarray(moving.ie, float)
    cont = np.asarray(static.ie_cnt, float)
    inner = (lam > lam[0] + 0.2) & (lam < lam[-1] - 0.2)
    mus = np.asarray(static.mus, float)
    assert i_static.shape[1] == mus.size == 5
    for ray, mu in enumerate(mus):
        spline = CubicSpline(lam, i_static[:, ray])
        err = np.abs(i_moving[inner, ray] - spline(lam[inner] * (1 + mu * vz / C_KMS)))
        flipped = np.abs(i_moving[inner, ray] - spline(lam[inner] * (1 - mu * vz / C_KMS)))
        assert np.max(err / cont[inner, ray]) < 1e-3, mu
        if mu > 0.5:
            assert np.max(flipped / cont[inner, ray]) > 0.1, mu


def test_constant_vz_leaves_the_departures_unchanged(tmp_path):
    """NLTE Mg on the frequency grid of atom.mg86g refined eightfold (nq_fac=8).

    Galilean invariance holds for the line transfer, but the continuum scattering term
    sigma (J - B) is evaluated in the observer frame, which breaks it at second order in
    vz in the outermost layers: 4e-4 at 1 km/s here. The default grid adds its own
    sampling error of a few percent at the top, which is why it is refined. Below
    log tau_500 = 0 the departures agree to the convergence limit.
    """
    nlte = f"""&atom_params         atom_file='{ATOM}' nq_fac=8 /
&nlte_params         convlim=1d-5 /
&rt_iteration_params maxiter=400 /
&rt_output_params    maxiter=0 solver='none' /
"""

    def departures(run):
        _, dep = run.get_depart()
        b = np.array([np.asarray(d, float).ravel() for d in dep]).T
        return b, np.log10(np.asarray(run.atmos.tau, float).ravel())

    b0, ltau = departures(_run(tmp_path, 'static', None, nlte))
    b1, _ = departures(_run(tmp_path, 'moving', 1.0, nlte))
    rel = np.abs(b1 / b0 - 1)
    assert rel.max() < 2e-3
    assert rel[:, ltau > 0].max() < 1e-5
