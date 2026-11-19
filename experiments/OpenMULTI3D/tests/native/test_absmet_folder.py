"""``&opac_table_params absmet_folder`` gives a run the table of the folder whose [Fe/H]
is nearest the run's FeH, the more metal-poor one on a tie, since a table carries no
metallicity scaling. dispatch.x records its choice as ``absmet_file`` in ``params.nml``.

Needs the built ``dispatch.x``, the solar MARCS model, ``atom.ba06`` and the absmet tables
of [Fe/H] = 0, -1 and -2.
"""
import os
import subprocess
from pathlib import Path

import pytest

import tumult

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_'
             'r+0.00_s+0.00.mod')
ATOM = Path('input/atoms/atom.ba06')
TABLES = {feh: Path(f'input/absmet_files/m{feh:+.2f}a{alpha:+.2f}c+0.00n+0.00o{alpha:+.2f}'
                    f'r+0.00s+0.00v+1.00.bin')
          for feh, alpha in ((0.0, 0.0), (-1.0, 0.4), (-2.0, 0.4))}

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file()
             and all((ROOT / p).is_file() for p in (MARCS, ATOM, *TABLES.values()))),
        reason='needs dispatch.x, the solar MARCS model, atom.ba06 and the absmet tables'),
]


@pytest.mark.parametrize('feh, table_feh', [(-0.4, 0.0), (-1.4, -1.0), (-1.5, -2.0)])
def test_the_table_nearest_the_run_feh(tmp_path, feh, table_feh):
    folder = tmp_path / 'absmet'
    folder.mkdir()
    for path in TABLES.values():
        (folder / path.name).symlink_to(ROOT / path)
    (folder / 'README.txt').write_text('not a table')
    (tmp_path / 'input').symlink_to(ROOT / 'input')
    (tmp_path / 'r.nml').write_text(f"""
&io_params datadir='{tmp_path / 'data'}' do_trace=F /
&atmos_params atmos_file='{MARCS}' FeH={feh} /
&atom_params atom_file='{ATOM}' /
&nlte_params convlim=1d-2 /
&rt_iteration_params maxiter=1 /
&opac_table_params absmet_folder='{folder}/' /
""")
    done = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=tmp_path, capture_output=True,
                          text=True, env=dict(os.environ, OMP_NUM_THREADS='2'))
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    nml = tumult.read(next((tmp_path / 'data').rglob('params.nml')).parent).nml
    assert nml['opac_table_params']['absmet_file'].strip() == str(folder / TABLES[table_feh].name)
