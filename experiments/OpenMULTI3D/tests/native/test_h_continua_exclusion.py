"""An H model atom replaces the background H continua unless exclude_h_continua is off.

HOP's opacity table carries the bound-free continua of hydrogen n = 1..5 as ``H BF``
entries. With ``exclude_trace_cont`` an H model atom with continua from these levels
takes them over, so that their opacity is counted once; with occupation probabilities
(the default) the table entries keep only their dissolved-level part redward of the
series limit. ``exclude_h_continua=F`` keeps the table continua fully active alongside
the atom's, as codes do that leave them in the background.

Needs the built ``dispatch.x``, the solar MARCS model and the public ``atom.h20``
(neither tracked in the repository).
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/public/atom.h20')
H_CONTINUA = {f'H BF G{n}' for n in range(1, 6)}

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM))),
        reason='needs dispatch.x, the p5777 MARCS model and atom.h20'),
]


def replaced_h_continua(tmp_path, exclude_h_continua):
    """The table H continua that a one-iteration h20 run excludes or reduces."""
    work = tmp_path / f'exclude_{exclude_h_continua}'
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    flag = 'T' if exclude_h_continua else 'F'
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&console_params level=2 /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' exclude_h_continua={flag} /
&nlte_params convlim=1d-3 /
&rt_iteration_params maxiter=1 /
&rt_output_params maxiter=0 /
""")
    env = dict(os.environ, OMP_NUM_THREADS='4')
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, env=env,
                            capture_output=True, text=True, timeout=900)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    pattern = r'(?:excluding transition|keeping only the dissolved-level part of) (H BF G\d)'
    return set(re.findall(pattern, result.stdout))


def test_h_atom_replaces_the_table_h_continua(tmp_path):
    assert replaced_h_continua(tmp_path, True) == H_CONTINUA


def test_exclude_h_continua_off_keeps_the_table_h_continua(tmp_path):
    assert replaced_h_continua(tmp_path, False) == set()
