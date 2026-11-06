"""A 1D NLTE run does not depend on the thread count.

The wavelength decomposition (``n_windows``) defaults to the thread count, so it is fixed
here; with it fixed, the populations and the spectrum are a property of the model
alone. The i_mu twins of each depth patch run each iteration concurrently and exchange
their hemisphere mean intensities, and the rt tasks of all frequency windows add into
the same profile normalisation. While the twin exchange depended on the order in
which tasks happened to run, repeated 8-thread runs differed by up to 1e-2 in the
populations of the upper patches and 4e-3 in the flux. The profile normalisation is a
real(4) sum in completion order, which leaves differences of about 1e-6; the tolerance
admits those and nothing near the size of a schedule-dependent iteration.

Needs the built ``dispatch.x`` and, from the input bundle (see ``README_INSTALL.md``),
the solar MARCS model, the ``mg86g`` model atom and ``SIU_masterline_Bergemann2026.txt``.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import tumult

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/atom.mg86g')
LINELIST = Path('input/linelists/SIU_masterline_Bergemann2026.txt')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM, LINELIST))),
        reason='needs dispatch.x, the solar MARCS model, atom.mg86g and SIU_masterline_Bergemann2026.txt'),
]


def _run(tmp_path, threads):
    """Five NLTE iterations and a three-iteration output scattering loop over 8 depth
    patches; returns the run's save directory."""
    work = tmp_path / f't{threads}'
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='{work / 'data'}' gb_step=6.0 do_trace=F /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' /
&nlte_params convlim=1d-2 decouple_continuum=T /
&rt_iteration_params maxiter=5 /
&rt_output_params maxiter=3 /
&run_params n_windows=8 /
&linelist_params linelist_file='{LINELIST}' line_scattering=T /
&spectrum_params daa=0.05 aa_blue=5160 aa_red=5192 /
""")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, env=env,
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    return work / 'data' / 'r'


def test_threaded_run_matches_serial_run(tmp_path):
    serial, threaded = (tumult.read(_run(tmp_path, n)) for n in (1, 8))
    for name in ('pop_Mg', 'flux'):
        if name == 'flux':
            a, b = (np.asarray(run.flux, float) for run in (serial, threaded))
        else:
            a, b = (np.asarray(run.read_patch_save(name, zdim=2), float)
                    for run in (serial, threaded))
        assert a.shape == b.shape and a.size > 0
        np.testing.assert_allclose(b, a, rtol=1e-4, err_msg=name)
