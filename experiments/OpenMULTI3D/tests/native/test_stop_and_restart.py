"""A stop request ends the NLTE iterations with every depth patch's departures written, and a
run that reads those departures (``read_departures``) and iterates on reaches the departures
of an uninterrupted run.

The request is ``<run_dir>/stop_nlte.flag``, written here once the first periodic save
(``save_every=1``) shows that the iterations are under way. The restart continues from the
saved b rather than from LTE, so it converges to the reference within ``convlim``; the
saved b are single precision.

Needs the built ``dispatch.x``, the solar MARCS model and ``atom.mg86g``.
"""
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

import tumult

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


def _namelist(work, atmos, extra=''):
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='{work / 'data'}' do_trace=F /
&atmos_params {atmos} /
&atom_params atom_file='{ATOM}' relative_abundance=T /
&composition_params abund_file='input/abund/abund_magg22' /
&opac_table_params absmet_file='input/absmet_files/absmet' /
&nlte_params convlim=1d-5 /
&rt_iteration_params maxiter=400 /
&spectrum_params aa_blue=5170 aa_red=5175 R_sampling=2e5 /
{extra}
""")


def _start(work):
    return subprocess.Popen([str(EXECUTABLE), 'r.nml'], cwd=work, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='2'))


def _b(run_dir):
    _, dep = tumult.read(run_dir).get_depart()
    return np.array([np.asarray(d, float).ravel() for d in dep])


def _rows(meta):
    return sum(1 for line in meta.read_text().splitlines() if line and not line.startswith('*'))


def test_stopped_run_restarts_to_the_reference(tmp_path):
    ref = tmp_path / 'ref'
    ref.mkdir()
    _namelist(ref, f"atmos_format='marcs' atmos_file='{MARCS}'")
    log = _start(ref).communicate(timeout=300)[0]
    assert 'OpenMULTI3D terminated successfully' in log, log[-2000:]

    stop = tmp_path / 'stop'
    stop.mkdir()
    _namelist(stop, f"atmos_format='marcs' atmos_file='{MARCS}'", '&output_params save_every=1 /')
    run_dir = stop / 'data' / 'r'
    proc = _start(stop)
    deadline = time.monotonic() + 120
    while not (run_dir / 'save' / 'pop_Mg_meta.txt').exists():
        assert proc.poll() is None and time.monotonic() < deadline, 'no periodic save appeared'
        time.sleep(0.02)
    (run_dir / 'stop_nlte.flag').write_text('')
    log = proc.communicate(timeout=300)[0]
    assert proc.returncode == 0 and 'OpenMULTI3D stopped on request' in log, log[-2000:]
    assert _rows(run_dir / 'save' / 'pop_Mg_meta.txt') == _rows(run_dir / 'save' / 'atmos_meta.txt')

    restart = tmp_path / 'restart'
    restart.mkdir()
    _namelist(restart, f"atmos_format='save' atmos_file='{run_dir}' read_departures=T")
    log = _start(restart).communicate(timeout=300)[0]
    assert 'OpenMULTI3D terminated successfully' in log, log[-2000:]

    b_ref = _b(ref / 'data' / 'r')
    assert np.max(np.abs(_b(run_dir) / b_ref - 1)) > 1e-3          # stopped well short
    assert np.max(np.abs(_b(restart / 'data' / 'r') / b_ref - 1)) < 1e-4
