"""A run named by an absolute path is written to that path.

``tm.launch`` reports ``<cwd>/<datadir>/<name>`` as the run folder, which for an absolute
name is the name itself; dispatch.x has to write the run there too, so that a synthesis
cache given as an absolute path finds its runs. The run is the 1D solar Ba demo of
``examples/sun_ba_1d.nml`` (well under a second).

Needs the built ``dispatch.x``.
"""
from pathlib import Path

import pytest

import tumult
from tumult.launch import run_complete

ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.skipif(not (ROOT / 'dispatch.x').is_file(), reason='needs dispatch.x')


def test_an_absolute_run_name_is_the_run_folder(tmp_path):
    work = tmp_path / 'work'
    work.mkdir()
    for name in ('dispatch.x', 'input'):
        (work / name).symlink_to(ROOT / name)
    name = tmp_path / 'elsewhere' / 'ba'
    folder, ran = tumult.launch(str(name), default=str(ROOT / 'examples' / 'sun_ba_1d.nml'),
                                cwd=work, verbose=0, threads=1)
    assert ran and folder == name
    assert run_complete(name)
    assert not (work / 'data').exists() or not any((work / 'data').rglob('abund.txt'))
