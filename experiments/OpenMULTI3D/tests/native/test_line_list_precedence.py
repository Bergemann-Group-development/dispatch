"""With a line list and an explicit output grid, the final spectrum takes the lines of the
atom's species from the line list, computed with the atom's departure coefficients.

The atom's own transitions drive the iterations only. The check is that the NLTE shift of
the Ba II 4554 core, relative to the LTE spectrum of the same run, is the one the atom's
own line gives (``output_from_line_list=F``): the list line, which has no level indices, is
matched to the atom's levels (``match_nlte_levels``) and takes its departures. Without that
flag only lines with level indices are in NLTE, so the list line stays in LTE. The wings
differ, since the list carries hyperfine and isotope components and its own damping.

Needs the built ``dispatch.x``, the solar MARCS model, ``atom.ba06`` and ``master_line.txt``.
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
ATOM = Path('input/atoms/atom.ba06')
LINELIST = Path('input/linelists/master_line.txt')

pytestmark = [
    pytest.mark.requires_input_data,
    pytest.mark.skipif(
        not (EXECUTABLE.is_file() and all((ROOT / p).is_file() for p in (MARCS, ATOM, LINELIST))),
        reason='needs dispatch.x, the solar MARCS model, atom.ba06 and master_line.txt'),
]


def _run(tmp_path, name, extra, atom=True, rt='', ll=''):
    work = tmp_path / name
    work.mkdir()
    (work / 'input').symlink_to(ROOT / 'input')
    (work / 'r.nml').write_text(f"""
&io_params datadir='{work / 'data'}' gb_step=6.0 do_trace=F /
&atmos_params atmos_format='marcs' atmos_file='{MARCS}' /
{f"&atom_params atom_file='{ATOM}' {extra} /" if atom else ''}
{"&nlte_params convlim=1d-2 /" if atom else ''}
&rt_output_params {rt} /
&linelist_params linelist_file='{LINELIST}' {ll} /
&spectrum_params R_sampling=5e5 aa_blue=4549 aa_red=4559 /
""")
    result = subprocess.run([str(EXECUTABLE), 'r.nml'], cwd=work, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS='2'), timeout=300)
    assert result.returncode == 0 and 'OpenMULTI3D terminated' in result.stdout, result.stdout[-2000:]
    return result.stdout, tumult.read(work / 'data' / 'r')


def _core_shift(run):
    return float(np.min(run.crop().ydata) - np.min(run.crop(LTE=True).ydata))


def test_list_line_takes_the_atoms_departures(tmp_path):
    log, from_list = _run(tmp_path, 'list', '', ll='match_nlte_levels=T')
    old_log, from_atom = _run(tmp_path, 'atom', 'output_from_line_list=F')
    assert 'come from the line list' in log
    assert 'excluding all Ba lines from line list' in old_log
    assert _core_shift(from_list) == pytest.approx(_core_shift(from_atom), abs=5e-3)
    assert _core_shift(from_list) < -0.02


def test_list_line_without_level_indices_stays_in_lte(tmp_path):
    _, unmatched = _run(tmp_path, 'unmatched', '')
    assert _core_shift(unmatched) == 0


def _equivalent_width(spectrum):
    y, x = (np.asarray(getattr(a, 'value', a)) for a in (spectrum.ydata, spectrum.xdata))
    return float(np.trapezoid(1 - y, x))


def test_lte_set_with_sampled_line_list_matches_plain_lte(tmp_path):
    """The LTE set of an NLTE run keeps the atom's species in the sampled-table path.

    ``sample_linelist`` leaves the departure species out of the sampled table; the LTE
    set must add them back with b=1, so that it reproduces a plain LTE run on the same
    line list and wavelengths.
    """
    _, plain = _run(tmp_path, 'plain', '', atom=False)
    _, sampled = _run(tmp_path, 'sampled', '', rt='sample_linelist=T sample_back_opac=T')
    reference = _equivalent_width(plain.crop())
    assert reference > 0
    assert _equivalent_width(sampled.crop(LTE=True)) == pytest.approx(reference, rel=2e-2)
