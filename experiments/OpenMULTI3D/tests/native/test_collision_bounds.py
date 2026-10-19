"""CH0 warnings use the physical rate coefficient for either ch0_hmip setting.

The Ba demo supplies the radiative atom. Its collision block is replaced with one
CH0 record, so the warning count and offending transition are unambiguous. The
valid example q proportional to T^(-1/2) is the Maxwellian average of sigma(E)
proportional to 1/E, with zero threshold. No atomic level energy enters the bound.
"""
import math
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'

pytestmark = pytest.mark.skipif(not EXECUTABLE.is_file(), reason='needs dispatch.x')


@pytest.mark.parametrize('hmip', [False, True])
@pytest.mark.parametrize('case', ['physical', 'falling', 'within_tolerance', 'unconverted'])
def test_ch0_maxwellian_warning(tmp_path, hmip, case):
    temps = (3000.0, 6000.0)
    # H- electron affinity divided by k_B, in kelvin.
    affinity_k = 0.754 / 8.617333262e-5
    if case in ('physical', 'unconverted'):
        rates = [1e-12 * (t / temps[0])**-0.5 for t in temps]
    else:
        rates = [1e-12 * (t / temps[0])**-1.5 for t in temps]
        rates[1] /= 1.20 if case == 'falling' else 1.02
    # The unconverted table retains the affinity factor even with the option off.
    if hmip or case == 'unconverted':
        rates = [q * math.exp(affinity_k / t) for q, t in zip(rates, temps)]

    (tmp_path / 'input').symlink_to(ROOT / 'input')
    atom = (ROOT / 'input/atoms/atom.ba06').read_text().split('GENCOL')[0]
    atom += ('GENCOL\nTEMP\n2 3000 6000\nCH0\n2 6 '
             + ' '.join(f'{q:.12e}' for q in rates) + '\nEND\n')
    (tmp_path / 'atom.ba_test').write_text(atom)
    (tmp_path / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&console_params level=2 /
&run_params mem_limit_gb=1 /
&atmos_params atmos_file='input/atmos/sun.mod' /
&atom_params atom_file='atom.ba_test' ch0_hmip={'.true.' if hmip else '.false.'} /
&nlte_params convlim=1d-2 /
&rt_iteration_params maxiter=1 /
&rt_output_params maxiter=1 /
""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=tmp_path,
                            env=dict(os.environ, OMP_NUM_THREADS='1'),
                            capture_output=True, text=True, timeout=120)
    output = result.stdout + result.stderr
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in output, output[-4000:]
    warned = case == 'falling' or (case == 'unconverted' and not hmip)
    assert ('collision record CH0 2 6 falls faster' in output) == warned, output
    assert ('1 of 1 collision records fall faster' in output) == warned, output
