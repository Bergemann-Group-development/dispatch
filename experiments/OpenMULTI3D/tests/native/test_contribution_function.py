"""The depression contribution function integrates to the computed line depression.

Along a ray, the difference of the transfer equations of the spectrum and of the
continuum gives D(0) = I_c(0) - I(0) = int (chi_l I_c - eta_l) exp(-tau_nu) ds exactly,
so ``Run.get_cntrbf`` integrated over log tau_500 has to reproduce F_cnt - F for the
flux and I_cnt - I for each ray, whatever the line strength. Ba II 4554 is a strong
line whose core still has a large opacity at the top of the model, where the optical
depth along the ray has to start from zero. The integral over log tau_500 uses the
PCHIP quadrature of the tau_500 scale itself (``hopacity.pchip_integral``); on the
64-point mesh it is within 0.75 % of the depression for the flux and both rays, the
core at mu = 0.2 included. A missing 1/mu, wrong flux weights or per-Hz units fail the
1.5 % tolerance by far.

Needs the built ``dispatch.x`` and the tracked solar MARCS model, Ba atom and
``master_line.txt``.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import tumult
from hopacity import pchip_integral

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/atom.ba06')
LINELIST = Path('input/linelists/master_line.txt')

pytestmark = pytest.mark.skipif(
    not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM, LINELIST))),
    reason='needs dispatch.x, the solar MARCS model, atom.ba06 and master_line.txt')


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    work = tmp_path_factory.mktemp('cntrbf')
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='{work / 'data'}' gb_step=6.0 do_trace=F /
&atmos_params atmos_file='{MARCS}' nz=64 /
&atom_params atom_file='{ATOM}' /
&nlte_params convlim=1d-2 /
&linelist_params linelist_file='{LINELIST}' /
&spectrum_params R_sampling=5e5 aa_blue=4552.5 aa_red=4555.5 /
&output_params save_cntrbf=T rt_to_file=F /
""")
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='256M')
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, env=env,
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    return tumult.read(work / 'data' / 'r')


@pytest.mark.parametrize('imu, rtol', [(None, 0.015), (0, 0.015), (-1, 0.015)])
def test_integral_reproduces_depression(run, imu, rtol):
    ltau, cf = run.get_cntrbf(imu=imu)
    integral = np.abs(pchip_integral(cf, ltau, axis=-1)[..., -1])
    if imu is None:
        depression = run.flux_cnt - run.flux
    else:
        depression = np.mean(run.ie_cnt[:, run.imus == run.imus[imu]]
                             - run.ie[:, run.imus == run.imus[imu]], axis=1)
    line = depression > 0.02 * depression.max()
    assert line.sum() > 10
    np.testing.assert_allclose(integral[line], depression[line], rtol=rtol)
