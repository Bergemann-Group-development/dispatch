"""Each window of a window file gives the spectrum it gives alone, however the output
patches split it.

The output windows share ``n_windows`` output patches in proportion to their size, with
no patch below 16 points unless its window is smaller (``set_n_windows_out``). The run uses
the sampled opacity tables in 1D, where a window's lines depend on two things: every
rt_patch of an output window must attach the shared line-list table, also when the stored
1D background lets it return early, and no patch may hold a table of one wavelength, which
the interpolation cannot bracket.

Needs the built ``dispatch.x``, the solar MARCS model and ``master_line.txt``.
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
LINELIST = Path('input/linelists/master_line.txt')
# Mg I b2, Mg I 5528 and Mg I 5711 with blends; 117, 110 and 72 points at R = 200000
WINDOWS = [(5171.184, 5174.184), (5526.905, 5529.905), (5710.088, 5712.088)]

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, LINELIST))),
        reason='needs dispatch.x, the solar MARCS model and master_line.txt'),
]


def _run(tmp_path, n_windows):
    work = tmp_path / f'nw{n_windows}'
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'w.txt').write_text(''.join(f'{a} {b}\n' for a, b in WINDOWS))
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&atmos_params atmos_file='{MARCS}' /
&run_params n_windows={n_windows} /
&rt_iteration_params maxiter=0 /
&rt_output_params maxiter=1 sample_back_opac=T sample_linelist=T /
&linelist_params linelist_file='{LINELIST}' /
&spectrum_params window_file='w.txt' R_sampling=200000.0 /
""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='4'), timeout=300)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    run = tumult.read(work / 'data' / 'r')
    sizes = [int(row.split()[2]) for row in
             (run.sfolder / 'flux_meta.txt').read_text().splitlines()[1:]]
    return run, sizes


def test_windows_do_not_depend_on_the_output_split(tmp_path):
    whole, whole_sizes = _run(tmp_path, 1)
    split, split_sizes = _run(tmp_path, 64)
    assert len(whole_sizes) == len(WINDOWS)
    assert len(split_sizes) > 2 * len(WINDOWS)
    assert min(split_sizes) >= 16
    for seg in whole.segments:
        window = (seg['blue'], seg['red'])
        ref = whole.crop(window=window).ydata
        flux = split.crop(window=window).ydata
        assert len(flux) == len(ref) > 50
        assert np.all(np.isfinite(flux))
        assert np.min(ref) < 0.5          # the window's lines are there
        np.testing.assert_allclose(flux, ref, rtol=1e-4)
