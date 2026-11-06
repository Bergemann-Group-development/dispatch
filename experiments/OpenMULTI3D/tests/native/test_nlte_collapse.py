"""The public Na atom neither collapses nor stalls in a cool 1D atmosphere.

``na_qmh`` has 149 Na II levels that lie 33 to 47 eV above the Na II ground and a Na III
level at the top.  In a cool atmosphere their LTE populations are 1e-60 to 1e-80 of the
Na, so two single-precision limits are reached.  The recombination coefficient into the
Na III level is the clamped Boltzmann factor (``tiny(1.0)``) times a LTE ratio of 1e65 or
more, which overflows ``real(4)``.  The population floor ``min_pop_frac`` then lifts the
Na II levels to a departure coefficient of 1e28 to 1e50, which the single-precision
opacities turn into ``Inf * 0``.  Either one made the statistical-equilibrium solution
``NaN``, the floor replaced it, and every level of the cell sat at ``1e-30`` of the atom,
with the run reporting a zero change and ``terminated successfully``.

A converged run satisfies particle conservation per cell (the populations sum to the LTE
total) and keeps the Na I ground state near LTE.  The collapsed solution sums to 1e-34 of
the total.  The 12 iterations suffice: the collapse had happened by iteration 3.

Needs the built ``dispatch.x``, the cool dwarf MARCS model, ``SIU_masterline_Bergemann2026.txt`` and the
public ``atom.na_qmh`` (neither tracked in the repository).
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from tumult import atom_rates

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p4000_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/public/atom.na_qmh')
LINELIST = Path('input/linelists/SIU_masterline_Bergemann2026.txt')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM, LINELIST))),
        reason='needs dispatch.x, the p4000 MARCS model, atom.na_qmh and SIU_masterline_Bergemann2026.txt'),
]


def _run(tmp_path):
    work = tmp_path / 'run'
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' abundance=6.33 /
&nlte_params convlim=1d-5 /
&rt_iteration_params maxiter=12 /
&linelist_params linelist_file='{LINELIST}' only_elements='Na' /
&spectrum_params daa=0.05 aa_blue=5682 aa_red=5690 /
&output_params save_net_rates=T /
""")
    env = dict(os.environ, OMP_NUM_THREADS='4')
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, env=env,
                            capture_output=True, text=True, timeout=900)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    return work / 'data' / 'r'


def test_cool_atmosphere_keeps_its_particle_budget(tmp_path):
    rates = atom_rates.read_rates(_run(tmp_path), 'Na')
    pop, nstar = rates['pop'], rates['nstar']
    np.testing.assert_allclose(pop.sum(0) / nstar.sum(0), 1.0, rtol=1e-3)
    ground_b = pop[0] / nstar[0]
    assert ground_b.min() > 0.1 and ground_b.max() < 10, ground_b
