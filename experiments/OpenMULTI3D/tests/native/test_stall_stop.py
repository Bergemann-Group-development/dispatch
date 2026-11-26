"""NLTE iterations that stop improving end early: with a convlim below the precision floor, the
largest change stops reaching new minima, and after ``stall_iter`` such iterations the run ends
as after a stop request, with departures written and the output pass run.

Needs the built ``dispatch.x``, the solar MARCS model and ``atom.mg86g``.
"""
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/atom.mg86g')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM))),
        reason='needs dispatch.x, the solar MARCS model and atom.mg86g'),
]


def test_stalled_iterations_stop_before_maxiter(tmp_path):
    (tmp_path / 'input').symlink_to(ROOT / 'input')
    (tmp_path / 'r.nml').write_text(f"""
&io_params datadir='{tmp_path / 'data'}' do_trace=F /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' relative_abundance=T /
&composition_params abund_file='input/abund/abund_magg22.txt' /
&opac_table_params absmet_file='input/absmet_files/m+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00v+1.00.bin' /
&nlte_params convlim=1d-12 convmax=0 stall_iter=10 /
&rt_iteration_params maxiter=200 /
&spectrum_params aa_blue=5170 aa_red=5175 R_sampling=2e5 /
""")
    proc = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=tmp_path, capture_output=True, text=True,
                          timeout=300, env=dict(os.environ, OMP_NUM_THREADS='2'))
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0 and 'NLTE iterations did not converge: they stalled after' in log, log[-2000:]
    save = tmp_path / 'data' / 'r' / 'save'
    assert (save / 'pop_Mg_meta.txt').is_file()
    assert any(save.glob('flux_*.bin'))
