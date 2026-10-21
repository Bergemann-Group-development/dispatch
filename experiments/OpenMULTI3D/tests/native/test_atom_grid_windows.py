"""A window file without a sampling grid crops the model atom's own output grid.

Without ``R_sampling``, ``daa``, ``nnu`` or ``lam_file`` the output pass runs on the
wavelengths of the atom's transitions, and the windows of a ``window_file`` keep only the
points inside them. The iteration pass is the same either way, so every kept point has the
flux of the uncropped run, to the 1e-6 by which two threaded NLTE solves differ in the
single-precision flux.

Needs the built ``dispatch.x`` and the solar MARCS model.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import tumult

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'
MARCS = Path('input/atmos/sun.mod')
ATOM = Path('input/atoms/atom.ba06')
# Ba II 4554 and 5853 (air)
WINDOWS = [(4553.0, 4555.0), (5852.7, 5854.7)]

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM))),
        reason='needs dispatch.x, the solar MARCS model and atom.ba06'),
]


def _run(tmp_path, name, spectrum):
    work = tmp_path / name
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'w.txt').write_text(''.join(f'{a} {b}\n' for a, b in WINDOWS))
    (work / 'r.nml').write_text(f"""
&io_params datadir='data' do_trace=F /
&atmos_params atmos_file='{MARCS}' /
&atom_params atom_file='{ATOM}' /
&nlte_params convlim=1e-3 /
&rt_output_params maxiter=0 /
{spectrum}
""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='4'), timeout=300)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    return tumult.read(work / 'data' / 'r')


def test_window_file_crops_the_atom_grid(tmp_path):
    whole = _run(tmp_path, 'whole', '')
    cropped = _run(tmp_path, 'cropped', "&spectrum_params window_file='w.txt' /")
    lam = np.asarray(cropped.crop(window=(4000.0, 7000.0)).xdata.value)
    inside = np.zeros(lam.size, bool)
    for a, b in WINDOWS:
        inside |= (lam >= a - 0.01) & (lam <= b + 0.01)
    assert inside.all()
    for window in WINDOWS:
        ref, out = whole.crop(window=window), cropped.crop(window=window)
        assert out.xdata.size == ref.xdata.size > 10
        np.testing.assert_array_equal(out.xdata.value, ref.xdata.value)
        assert np.min(ref.ydata) < 0.9          # the window's line is there
        np.testing.assert_allclose(out.ydata, ref.ydata, rtol=1e-5)
