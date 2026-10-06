"""The lines read for a window file depend only on the union of its windows.

A window file of several narrow windows within ``dlam`` of each other merges to one line-list
interval (``merge_intervals``), and the Intel builds, which compile with
``-assume norealloc_lhs``, once read no line at all from that case because the merged array
kept the size of the unmerged one. The checks compare the number of lines read with the
count of the equivalent single windows. They run on whatever compiler built ``dispatch.x``;
the regression only shows on an ifort build or a gfortran one with ``-fno-realloc-lhs``.

Needs the built ``dispatch.x``, the solar MARCS model and ``master_line.txt``.
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
LINELIST = Path('input/linelists/master_line.txt')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, LINELIST))),
        reason='needs dispatch.x, the solar MARCS model and master_line.txt'),
]


def _lines_read(tmp_path, windows):
    work = tmp_path / ('w' + '_'.join(f'{a:g}-{b:g}' for a, b in windows))
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'w.txt').write_text(''.join(f'{a} {b}\n' for a, b in windows))
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&atmos_params atmos_file='{MARCS}' /
&rt_iteration_params maxiter=0 /
&rt_output_params maxiter=3 /
&linelist_params linelist_file='{LINELIST}' /
&spectrum_params window_file='w.txt' R_sampling=200000.0 /
""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=600)
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]
    return int(re.search(r'^(\d+) of \d+ lines read from', result.stdout, re.M).group(1))


def test_windows_that_merge_read_the_lines_of_their_union(tmp_path):
    narrow = [(6150, 6166), (6292, 6308), (6357, 6370)]
    n = _lines_read(tmp_path, narrow)
    assert n > 0
    assert n == _lines_read(tmp_path, [(6150, 6370)])


def test_windows_that_stay_apart_read_the_lines_of_each(tmp_path):
    a, b = (4000, 4010), (7000, 7010)
    n = _lines_read(tmp_path, [a, b])
    assert n > 0
    assert n == _lines_read(tmp_path, [a]) + _lines_read(tmp_path, [b])
