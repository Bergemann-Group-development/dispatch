"""LTE opacities at individual (T, rho) points over the domain of the campaign atmospheres.

    python opacity_points.py points                         # write opac/points.npz
    python opacity_points.py shard FEH K N                  # shard K of N at [Fe/H]=FEH

``points`` samples (log T, log rho) uniformly (scrambled Sobol) inside the region the
atmospheres occupy: every depth point of the batch-1 MARCS models and every cell of the
first snapshot of each pilot Stagger model, accepted within RADIUS of an occupied point
in (log T / 0.05, log rho / 0.25) units, the spacing of HOP's own EOS tables.

``shard`` solves the EOS on its slice of the points at one composition (magg22 scaled to
[Fe/H], alpha enhanced as in the MARCS standard composition), then writes the total
extinction chi, its continuum part and the scattering part [1/cm, float32] on the HRS
grid of batch 1 (R=2e5), with the FGK line list of make_nml.py, vmic = 1 km/s.  Each
shard is one .npz under runs/opac/feh<FEH>/.

The line list is read through the windowed reader: the range reader crashes in the
frozen ifort build bin/5eb036cbb (build_linelist's ``blk_lo = [lo]`` under -assume
norealloc_lhs, fixed in a9b67c3e2).
"""
import sys
from pathlib import Path

import numpy as np

CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')
WORK = CAMPAIGN / 'work'
OUT = CAMPAIGN / 'runs/opac'
JOBS = Path(__file__).resolve().parent
N_POINTS = 8192          # per composition
RADIUS = 1.0             # acceptance radius in (dlogT/0.05, dlogrho/0.25) units
PILOT_SNAPSHOTS = ['t5777/t5777g44m00', 't60/t60g40m00', 't60/t60g40m20', 't55/t55g45m00_new',
                   't50/t50g25m00', 't45/t45g25m00_new', 't45/t45g25m20_new']


def load_make_nml():
    import importlib.util
    spec = importlib.util.spec_from_file_location('make_nml', JOBS / 'make_nml.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def marcs_t_rho(path):
    lines = Path(path).read_text().splitlines()
    i = next(j for j, line in enumerate(lines) if 'Number of depth points' in line)
    nz = int(lines[i].split()[0])
    temp = [float(lines[i + 3 + k].split()[4]) for k in range(nz)]
    rho = [float(lines[i + 4 + nz + k].split()[3]) for k in range(nz)]
    return np.log10(temp), np.log10(rho)


def stagger_t_rho(folder):
    atm = sorted((WORK / 'stagger2' / folder).glob('atm3d.*'))[0]
    n = 80 * 80 * 240
    data = np.fromfile(atm, dtype='<f4', count=6 * n).reshape(6, n)
    return np.log10(data[1]), np.log10(data[5])          # temp, rho


def occupied():
    mk = load_make_nml()
    cols = [marcs_t_rho(WORK / 'marcs' / p.name) for p, _ in mk.domain_models(WORK / 'marcs')]
    cols += [stagger_t_rho(f) for f in PILOT_SNAPSHOTS]
    lt = np.concatenate([c[0] for c in cols])
    lr = np.concatenate([c[1] for c in cols])
    return np.column_stack([lt / 0.05, lr / 0.25])


def make_points():
    from scipy.spatial import cKDTree
    from scipy.stats import qmc
    occ = occupied()
    tree = cKDTree(np.unique(np.round(occ, 2), axis=0))
    lo, hi = occ.min(axis=0) - RADIUS, occ.max(axis=0) + RADIUS
    sobol = qmc.Sobol(2, scramble=True, seed=20261002)
    kept = []
    while sum(len(k) for k in kept) < N_POINTS:
        u = qmc.scale(sobol.random(4096), lo, hi)
        d, _ = tree.query(u, distance_upper_bound=RADIUS)
        kept.append(u[np.isfinite(d)])
    pts = np.concatenate(kept)[:N_POINTS]
    temp, rho = 10 ** (pts[:, 0] * 0.05), 10 ** (pts[:, 1] * 0.25)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / 'points.npz', temp=temp, rho=rho,
             box_logT=np.array([lo[0], hi[0]]) * 0.05, box_logrho=np.array([lo[1], hi[1]]) * 0.25)
    print(f'{len(temp)} points: T {temp.min():.0f}-{temp.max():.0f} K, '
          f'rho {rho.min():.2e}-{rho.max():.2e} g/cm3 from {len(occ)} occupied states')


def alpha_for(feh):
    """[alpha/Fe] of the MARCS standard composition."""
    return float(np.clip(-0.4 * feh, 0.0, 0.4))


def shard(feh, k, n):
    import time
    from tumult.hop import HopEOS, wavelength_grid
    from tumult.physics import air2vac
    mk = load_make_nml()
    lists = [str(WORK / mk.MASTER / f) for f in mk.FGK_LINELISTS]
    rows = [tuple(float(x) for x in line.split())
            for line in (JOBS / 'windows/hrs.txt').read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    lam = np.concatenate([wavelength_grid(b, r, R_sampling=2e5) for b, r in rows])
    win = [(air2vac(b), air2vac(r)) for b, r in rows]
    dest = OUT / f'feh{feh:+.2f}'
    target = dest / f'shard_{k:04d}_of_{n:04d}.npz'
    if target.exists():
        return
    pts = np.load(OUT / 'points.npz')
    sel = slice(k, None, n)
    temp, rho = pts['temp'][sel], pts['rho'][sel]
    t0 = time.perf_counter()
    with HopEOS(feh=feh, alpha_enhancement=alpha_for(feh)) as hop:
        hop.load_state(temp, rho=rho, vmic=1.0)
        ne, pg = hop.ne, hop.pg
        chi_cont, _ = hop.opacity(lam)
        chi, scat = hop.opacity(lam, linelist=lists, dlam=5.0, windows=win)
    dest.mkdir(parents=True, exist_ok=True)
    np.savez(target, index=np.arange(len(pts['temp']))[sel],
             temp=temp, rho=rho, ne=ne, pg=pg, lam_vac=lam, chi=chi, chi_cont=chi_cont,
             scat=scat, feh=feh, alpha=alpha_for(feh), vmic=1.0)
    print(f'feh {feh:+.2f} shard {k}/{n}: {len(temp)} points in '
          f'{time.perf_counter() - t0:.1f} s')


if __name__ == '__main__':
    if sys.argv[1] == 'points':
        make_points()
    elif sys.argv[1] == 'shard':
        shard(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    else:
        sys.exit(__doc__)
