"""Standalone line-list pruning: the Python driver of hop_prune.x and its namelist.

The driver tests use a stand-in executable, so they run anywhere with a POSIX shell.
The tests of the executable itself are in ``tests/native/test_prune_executable.py``.
"""
import os
import stat
import threading
import time

import pytest

from tumult import prune
from tumult.synthesis.cancel import CancelledError, CancelToken

R = 300000.0
DLAM = 2.0


# --- the namelist ---------------------------------------------------------------------------
def test_prune_namelist_carries_the_runs_settings_and_the_node_overrides(tmp_path):
    base = {'atmos_params': {'atmos_file': 'other.mod', 'dims': 8, 'use_ne': False},
            'composition_params': {'abund_file': 'input/abund/abund_magg', 'abundances': 'Fe=7.0'},
            'opac_table_params': {'opac_table_vmic': 1.5},
            # a screen's element-only arm: its switch is that variant's, not the node's
            'linelist_params': {'linelist_folder': 'lists/', 'tolerance': 1e-3, 'dlam': 100.0,
                                'only_elements': 'Fe'},
            'spectrum_params': {'R': 1e5, 'aa_blue': 5000.0, 'in_air': False}}
    nml = prune.prune_namelist(tmp_path / 'o.txt', 'a.mod', ['m.txt'], tmp_path / 'w.txt',
                               R=R, dlam=DLAM, vmic=1.0, feh=-0.5, abundances='Fe=+0.4',
                               relative_abundances=True, base=base, use_abnd=True)
    assert nml['atmos_params'] == {'atmos_file': 'a.mod', 'use_ne': False, 'feh': -0.5,
                                   'vmic': 1.0, 'use_abnd': True}
    # the caller's rung replaces the base abundance override; the other settings survive
    assert nml['composition_params'] == {'abund_file': 'input/abund/abund_magg',
                                         'abundances': 'Fe=+0.4', 'relative_abundances': True}
    assert nml['opac_table_params'] == {'opac_table_vmic': 1.5}
    assert nml['linelist_params'] == {'tolerance': 1e-3, 'linelist_file': ['m.txt'],
                                      'dlam': DLAM}
    assert nml['spectrum_params'] == {'R': R, 'window_file': str(tmp_path / 'w.txt'),
                                      'in_air': False}
    assert nml['prune_params'] == {'output_file': str(tmp_path / 'o.txt'),
                                   'sample_count': prune.SAMPLE_COUNT}
    # a ladder that reads one species alone says so itself
    nml = prune.prune_namelist(tmp_path / 'o.txt', 'a.mod', ['m.txt'], tmp_path / 'w.txt',
                               R=R, dlam=DLAM, base=base, only_elements='Fe')
    assert nml['linelist_params']['only_elements'] == 'Fe'


# --- the driver, on a stand-in executable ----------------------------------------------------
def _fake_executable(tmp_path, body):
    script = tmp_path / 'fake_prune.sh'
    script.write_text('#!/bin/sh\n' + body)
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def test_cancellation_kills_the_child_and_publishes_nothing(tmp_path, monkeypatch):
    out = tmp_path / 'node.txt'
    # writes its partial like the real executable, then waits to be killed
    script = _fake_executable(tmp_path, f'echo partial > "{out}.$$.partial"\nsleep 30\n'
                                        f'mv "{out}.$$.partial" "{out}"\n')
    monkeypatch.setenv('HOP_PRUNE', str(script))
    nml_path = tmp_path / 'node.nml'
    nml_path.write_text('&prune_params /\n')
    token = CancelToken()
    result = {}

    def run():
        try:
            prune.run_prune(nml_path, out, cwd=tmp_path, cancel_token=token)
        except BaseException as exc:      # noqa: BLE001 - the outcome is the assertion
            result['exc'] = exc

    thread = threading.Thread(target=run)
    thread.start()
    deadline = time.time() + 10
    while not list(tmp_path.glob('node.txt.*.partial')) and time.time() < deadline:
        time.sleep(0.02)
    token.cancel()
    thread.join(10)
    assert not thread.is_alive()
    assert isinstance(result.get('exc'), CancelledError)
    assert not out.exists()
    assert not list(tmp_path.glob('node.txt.*.partial'))


def test_a_failing_child_raises_with_its_diagnostic(tmp_path, monkeypatch):
    script = _fake_executable(tmp_path, 'echo "hop_prune: file not found: x.mod" >&2\nexit 1\n')
    monkeypatch.setenv('HOP_PRUNE', str(script))
    nml_path = tmp_path / 'node.nml'
    nml_path.write_text('&prune_params /\n')
    with pytest.raises(RuntimeError, match='file not found: x.mod'):
        prune.run_prune(nml_path, tmp_path / 'node.txt', cwd=tmp_path)


def test_prune_batch_skips_cached_outputs_and_runs_the_rest_concurrently(tmp_path, monkeypatch):
    # each child records when it ran; two uncached nodes must overlap in time
    script = _fake_executable(tmp_path, 'date +%s.%N > "$1.started"\nsleep 1\n'
                                        'echo lines > "$(basename "$1" .nml).txt"\n')
    monkeypatch.setenv('HOP_PRUNE', str(script))
    cached = tmp_path / 'a.txt'
    cached.write_text('cached')
    jobs = [(cached, {'prune_params': {}}), (tmp_path / 'b.txt', {'prune_params': {}}),
            (tmp_path / 'c.txt', {'prune_params': {}})]
    t0 = time.time()
    outputs = prune.prune_batch(jobs, cwd=tmp_path)
    elapsed = time.time() - t0
    assert outputs == [cached, tmp_path / 'b.txt', tmp_path / 'c.txt']
    assert cached.read_text() == 'cached' and not (tmp_path / 'a.nml').exists()
    assert (tmp_path / 'b.txt').is_file() and (tmp_path / 'c.txt').is_file()
    if (os.cpu_count() or 1) > 1:
        assert elapsed < 1.9, 'the two uncached nodes did not run concurrently'
