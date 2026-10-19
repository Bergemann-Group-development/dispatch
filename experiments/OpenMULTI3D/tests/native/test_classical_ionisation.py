"""Classical ionisation couples a level only to the upper level of its continuum.

The tabulated ionisation rates of the public model atoms all lead from a level to the
upper level of its continuum, and so do DETAIL's generated ones. The classical recipes
(``e_ionisation``, ``h_ionisation``) follow them: with the recipes on, a pair of levels
without tabulated rates is coupled by ionisation exactly when the upper level ends the
lower level's continuum. ``o41f`` has excited O II levels and a continuum that ends in
one of them (2p3 2D), so both kinds of continuum are covered.

Needs the built ``dispatch.x``, the solar MARCS model and the public ``atom.o41f``
(neither tracked in the repository).
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from tumult import atom_rates
from tumult.atom import Atom
from tumult.atom_audit import continuum_level

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/public/atom.o41f')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM))),
        reason='needs dispatch.x, the p5777 MARCS model and atom.o41f'),
]


def test_classical_ionisation_leads_to_the_continuum_upper_level(tmp_path):
    work = tmp_path / 'run'
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' e_excitation='regemorter'
    e_ionisation='seaton' h_excitation='drawin' h_ionisation='drawin' /
&nlte_params convlim=1d-3 /
&rt_iteration_params maxiter=1 /
&rt_output_params maxiter=0 /
&output_params save_net_rates=T /
""")
    env = dict(os.environ, OMP_NUM_THREADS='4')
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, env=env,
                            capture_output=True, text=True, timeout=900)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]

    atom = Atom.from_old(ROOT / ATOM)
    ion = atom.lvl.ion.to_numpy()
    upper = continuum_level(atom)
    assert np.any((upper >= 0) & (upper != np.argmax(ion == 2)))   # an excited upper level
    tabulated = {tuple(sorted(p)) for p in zip(atom.col.i, atom.col.j)}
    cup = atom_rates.read_rates(work / 'data' / 'r', 'O')['cup']
    nk = len(ion)
    pair = 0
    for j in range(1, nk):                      # j-major, i < j
        for i in range(j):
            if ion[j] == ion[i] + 1 and (i, j) not in tabulated:
                coupled = np.any(cup[pair] > 0)
                assert coupled == (j == upper[i]), (i, j, upper[i])
            pair += 1
