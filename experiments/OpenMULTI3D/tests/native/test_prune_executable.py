"""Standalone line-list pruning: the hop_prune.x executable.

These tests need the standalone build (``make -C HOP/standalone``), the checked-in solar
MARCS model and the master line list, and are skipped without them. The Python driver
and its namelist are tested on a stand-in executable in
``tests/execution/test_prune_driver.py``.
"""
from pathlib import Path

import numpy as np
import pytest

from tumult import prune

ROOT = Path(__file__).resolve().parents[2]
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
    lam = np.concatenate([wavelength_grid(lo, hi, R_sampling=R) for lo, hi in WINDOWS])
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
    assert record['spectrum_params']['r_sampling'] == R
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
