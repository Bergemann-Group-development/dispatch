"""Phase-0 timing of LTE opacities at individual (T, rho) points through tumult.hop.

Samples points over the (T, rho) domain of FGK photospheres, solves the EOS there and
computes chi and its scattering part on the 4MOST HRS grid (R=2e5, the 1D benchmark's
grid), with the continuum alone and with the FGK line list of make_nml.py. Reports
points per second per process and the output volume per point. Nothing is saved.

The line list is read through the windowed reader (``windows=``): the plain range
reader (hop_init_linelist) crashes in the frozen ifort build bin/5eb036cbb, whose
build_linelist assigns ``blk_lo = [lo]`` to an unallocated array under -assume
norealloc_lhs (fixed in a9b67c3e2). HopEOS.opacity(nprocs>1) takes that range path in
its workers, so the parallel variant is left out; points parallelise as independent
processes instead.
"""
import time
from pathlib import Path

import numpy as np

from tumult.hop import HopEOS, wavelength_grid
from tumult.physics import air2vac

JOBS = Path(__file__).resolve().parent
CLONE = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D')
MASTER = CLONE / 'input/linelists/master_linelists'


def fgk_linelists():
    import importlib.util
    spec = importlib.util.spec_from_file_location('make_nml', JOBS / 'make_nml.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return [str(MASTER / f) for f in mod.FGK_LINELISTS]


def hrs_windows():
    """The HRS windows as (blue, red) air Angstrom."""
    return [tuple(float(x) for x in line.split())
            for line in (JOBS / 'windows/hrs.txt').read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def hrs_grid(R=2e5):
    return np.concatenate([wavelength_grid(b, r, R_sampling=R) for b, r in hrs_windows()])


def sample_points(n, seed=1):
    """Uniform in log T over 3000-12000 K and log rho over 1e-11-1e-6 g/cm3."""
    rng = np.random.default_rng(seed)
    return 10 ** rng.uniform(np.log10(3000), np.log10(12000), n), \
        10 ** rng.uniform(-11, -6, n)


def timed(label, fn, n_points, nlam):
    t0 = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t0
    print(f'{label:42s} {dt:8.2f} s  {n_points / dt:9.2f} points/s  '
          f'{n_points * nlam / dt / 1e6:8.2f} Mpoint-lam/s', flush=True)
    return out


def main():
    lam = hrs_grid()
    lists = fgk_linelists()
    win = [(air2vac(b), air2vac(r)) for b, r in hrs_windows()]
    print(f'HRS grid: {lam.size} wavelengths, {lam[0]:.1f}-{lam[-1]:.1f} A (vacuum)')
    for n in (64, 512):
        temp, rho = sample_points(n)
        with HopEOS(feh=0.0) as hop:
            timed(f'EOS, n={n}', lambda: hop.load_state(temp, rho=rho, vmic=1.0), n, 1)
            timed(f'continuum, n={n}', lambda: hop.opacity(lam), n, lam.size)
            chi, scat = timed(f'continuum+lines (incl. list read), n={n}',
                              lambda: hop.opacity(lam, linelist=lists, dlam=5.0, windows=win), n, lam.size)
            timed(f'continuum+lines (list in memory), n={n}',
                  lambda: hop.opacity(lam, linelist=lists, dlam=5.0, windows=win), n, lam.size)
            print(f'  output chi+scat float32: {chi.nbytes * 2 / n / 1e6:.3f} MB per point; '
                  f'median chi at 5174 A: '
                  f'{np.median(chi[:, np.searchsorted(lam, 5174.0)]):.3e} 1/cm')


if __name__ == '__main__':
    main()
