"""A run's model atom, read back from its ``atoms/`` folder.

dispatch.x writes the conversion of a MULTI atom there (``atom.ba06`` becomes
``ba06.txt``), and :meth:`tumult.run.Run.get_line` finds a transition in it whose levels
select departure coefficients: Ba II 4554 Å is 6s 2S - 6p 2P3/2, levels 1 and 4 counted
from 0. The run is the 1D solar demo of ``examples/sun_ba_1d.nml`` (well under a second).

Needs the built ``dispatch.x``.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import tumult

ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / 'dispatch.x'

pytestmark = pytest.mark.skipif(not EXECUTABLE.is_file(), reason='needs dispatch.x')


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    work = tmp_path_factory.mktemp('demo')
    (work / 'input').symlink_to(ROOT / 'input')
    nml = (ROOT / 'examples' / 'sun_ba_1d.nml').read_text().replace("datadir='data'", f"datadir='{work / 'data'}'")
    (work / 'demo.nml').write_text(nml)
    subprocess.run([str(EXECUTABLE), 'demo.nml'], cwd=work, check=True, capture_output=True,
                   env=dict(os.environ, OMP_NUM_THREADS='1'))
    return tumult.read(work / 'data' / 'demo')


def test_the_run_keeps_the_converted_atom(run):
    assert [p.name for p in (run.folder / 'atoms').iterdir()] == ['ba06.txt']


def test_a_line_s_levels_select_their_departure_coefficients(run):
    line = run.get_line(4554.033)
    assert (line.lower.i, line.upper.i) == (1, 4)
    for lvl in (line.lower, line.upper):
        np.testing.assert_array_equal(run.get_depart(lvl)[1], run.get_depart(lvl.i, 'Ba')[1])
