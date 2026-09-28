"""Standalone line-list pruning: the hop_prune.x executable and its Python driver.

The executable tests need the standalone build (``make -C HOP/standalone``) and the
checked-in solar MARCS model plus the master line list, and are skipped without them. The
driver tests use a stand-in executable, so they run anywhere.
"""
import os
import stat
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from tumult import prune
from tumult.synthesis.cancel import CancelledError, CancelToken

ROOT = Path(__file__).resolve().parent.parent
MARCS = ROOT / 'input' / 'atmos' / ('p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_'
                                    'o+0.00_r+0.00_s+0.00.mod')
MASTER = ROOT / 'input' / 'linelists' / 'master_line.txt'
WINDOWS = [(5000.0, 5002.5), (5015.0, 5017.5)]      # disjoint, air Angstrom
R = 300000.0
DLAM = 2.0

needs_executable = pytest.mark.skipif(
    not ((ROOT / prune.EXECUTABLE).is_file() and MARCS.is_file() and MASTER.is_file()),
    reason='needs HOP/standalone/build/hop_prune.x, the solar MARCS model and master_line.txt')


def _records(path):
    """``{(species, lam0, evi, loggf): strength}`` of a block-format line list."""
    out, species = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("'") and not line.strip().strip("'").split("'")[0].isdigit():
            species = line.strip().strip("'")
            continue
        parts = line.split()
        if len(parts) >= 9:
            key = (species, round(float(parts[0]), 4), round(float(parts[1]), 4),
                   round(float(parts[2]), 3))
            out[key] = float(parts[8])
    return out


def _prune(tmp_path, name, **overrides):
    out = tmp_path / f'{name}.txt'
    window_file = tmp_path / 'windows.txt'
    window_file.write_text('\n'.join(f'{lo} {hi}' for lo, hi in WINDOWS) + '\n')
    kwargs = dict(R=R, dlam=DLAM, vmic=1.0, feh=0.0)
    kwargs.update(overrides)
    nml = prune.prune_namelist(out, MARCS, [MASTER], window_file, **kwargs)
    prune.prune_batch([(out, nml)], cwd=ROOT)
    return out


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


# --- the executable ---------------------------------------------------------------------------
@needs_executable
def test_executable_selects_what_the_in_process_crop_selects(tmp_path):
    """Parity with the in-process HopEOS selection on the same states.

    The Python crop solves the EOS from rho alone; the executable follows the run's
    use_rho/use_ne settings, so the comparison sets use_ne=False. Sampling eight depths
    against every depth must also give the same set: the selection is a maximum over the
    states, and the sampling spans the structure.
    """
    import tumult as tm
    from tumult.hop import wavelength_grid
    from tumult.physics import air2vac
    from tumult.linestrength import read_marcs

    sampled = _prune(tmp_path, 'sampled', base={'atmos_params': {'use_ne': False}})
    every = _prune(tmp_path, 'every', base={'atmos_params': {'use_ne': False}},
                   sample_count=0)
    lam = np.concatenate([wavelength_grid(lo, hi, R=R) for lo, hi in WINDOWS])
    vac = [(float(air2vac(lo)), float(air2vac(hi))) for lo, hi in WINDOWS]
    temp, rho = read_marcs(MARCS)
    with tm.HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho, vmic=1.0)
        hop.write_linelist(tmp_path / 'python.txt', lam, str(MASTER), dlam=DLAM, vmic=1.0,
                           windows=vac)
    a, b, c = _records(sampled), _records(every), _records(tmp_path / 'python.txt')
    assert len(a) > 300
    assert set(a) == set(b) == set(c)
    record = prune.provenance(sampled)
    result = record['prune_result']
    assert result['ntransitions'] == len(a) and record['prune_params']['sample_count'] == 8
    assert result['sampled_depths'][0] == 1 and result['sampled_depths'][-1] == result['nz_model']
    assert len(result['blue_vac']) == 2 and record['linelist_params']['dlam'] == DLAM
    assert record['spectrum_params']['r'] == R
    assert result['use_rho_used'] is True and result['use_ne_used'] is False
    # boundary transitions: the wing contributors outside the windows are kept, nothing
    # farther than dlam from a window is
    lam0 = np.array([k[1] for k in a])
    inside = np.zeros(lam0.shape, bool)
    for lo, hi in WINDOWS:
        inside |= (lam0 >= lo - DLAM) & (lam0 <= hi + DLAM)
    assert inside.all()
    assert ((lam0 < WINDOWS[0][0]) | (lam0 > WINDOWS[0][1])).any()


@needs_executable
def test_strongest_rung_keeps_every_weaker_rungs_lines_of_its_element(tmp_path):
    default = _records(_prune(tmp_path, 'default'))
    strong = _records(_prune(tmp_path, 'strong', abundances='Fe=+1.0',
                             relative_abundances=True))
    iron = lambda recs: {k for k in recs if k[0] in ('Fe I', 'Fe II')}
    assert iron(default) <= iron(strong)
    assert all(strong[k] > default[k] for k in iron(default))
    composition = prune.provenance(tmp_path / 'strong.txt')['composition_params']
    assert composition['abundances'] == 'Fe=+1.0' and composition['relative_abundances'] is True


@needs_executable
def test_atmosphere_abundances_and_absolute_overrides_are_applied_in_run_order(tmp_path):
    """MARCS abundances replace the standard set when use_abnd is set, and an explicit
    override is applied on top of them -- the composition_type%init order a run uses."""
    marcs = _records(_prune(tmp_path, 'marcs', use_abnd=True))
    default = _records(_prune(tmp_path, 'default'))
    assert prune.provenance(tmp_path / 'marcs.txt')['prune_result']['use_abnd_used'] is True
    assert set(marcs) != set(default)
    # an absolute Fe override on top of the MARCS set: the same iron lines whichever standard
    # set sits underneath, since the override replaces the element's value
    marcs_fe = _records(_prune(tmp_path, 'marcs_fe', use_abnd=True, abundances='Fe=8.0',
                               relative_abundances=False))
    default_fe = _records(_prune(tmp_path, 'default_fe', abundances='Fe=8.0',
                                 relative_abundances=False))
    iron = lambda recs: {k: v for k, v in recs.items() if k[0] in ('Fe I', 'Fe II')}
    assert set(iron(marcs_fe)) == set(iron(default_fe))


@needs_executable
def test_invalid_inputs_fail_without_publishing(tmp_path):
    out = tmp_path / 'bad.txt'
    window_file = tmp_path / 'windows.txt'
    window_file.write_text('5000 5002.5\n')
    nml = prune.prune_namelist(out, MARCS, [tmp_path / 'missing.txt'], window_file,
                               R=R, dlam=DLAM, vmic=1.0, feh=0.0)
    with pytest.raises(RuntimeError, match='file not found'):
        prune.prune_batch([(out, nml)], cwd=ROOT)
    assert not out.exists() and not list(tmp_path.glob('bad.txt.*'))
    nml = prune.prune_namelist(out, tmp_path / 'nope.mod', [MASTER], window_file,
                               R=R, dlam=DLAM, vmic=1.0, feh=0.0)
    with pytest.raises(RuntimeError, match='File not found'):
        prune.prune_batch([(out, nml)], cwd=ROOT)
    assert not out.exists()
