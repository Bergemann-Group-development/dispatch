#!/usr/bin/env python3
"""<3D> and 1D MARCS model atmospheres at the benchmark stars' own parameters.

Both grids are interpolated the same way, on the log tau_500 scale of the models:

- <3D>: the campaign's averages ``work/avg3d8/<model>.atmos`` (``tumult.avg3d``, tau500
  surfaces, <T^4>^(1/4), top variant A; reports/2026-10-06_phase7.md), whose depth column
  is log tau_500;
- MARCS: the standard models (plane-parallel 1 Msun t01 for log g >= 3.5, spherical
  1 Msun t02 below, [alpha/Fe] of the standard composition), on their lgTau5 column.

At each depth, T, ln P_e, ln rho and vmic are fitted by a weighted least-squares plane in
(Teff/1000 K, log g, [Fe/H]) over the ``K`` nodes nearest the star in the scaled distance
(500 K, 0.5 dex, 0.5 dex), with weights 1/(d^2 + 0.1). This is first-order (Taylor)
interpolation, which the irregular <3D> grid allows where trilinear corners do not
exist; for the regular MARCS grid it reduces to the usual linear interpolation. The depth
range is the one all used nodes share. A node closer than ``EXACT`` to the star is used
as it is (the Sun on t5777g44m00). ``validate`` reports the leave-one-out error of
the <3D> interpolant over the grid nodes, the scale of its error at a star.

The output is dispatch.x's 'text' format with log tau_500 as the depth column, read with
``atmos_format='text' use_ne=F amr='T'`` like the campaign averages.

    python atmospheres.py AVG3D_DIR MARCS_DIR OUT_DIR       # writes <star>_av3d.atmos, <star>_marcs.atmos
    python atmospheres.py --validate AVG3D_DIR
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SCALE = np.array([0.5, 0.5, 0.5])          # Teff [kK], log g, [Fe/H] per unit distance
K_NODES = 8
EXACT = 0.1                                # closer than this, a node is used as it is
NDEPTH = 101
AV3D_NAME = re.compile(r't(?P<t>\d+)g(?P<g>\d\d)(?P<s>[mp])(?P<z>\d\d)$')
MARCS_NAME = re.compile(r'(?P<geo>[ps])(?P<t>\d{4})_g(?P<g>[+-]\d\.\d)_m(?P<m>[\d.]+)_t(?P<vt>\d\d)_st_'
                        r'z(?P<z>[+-]\d\.\d\d)_a(?P<a>[+-]\d\.\d\d)_c\+0\.00_n\+0\.00_o(?P<o>[+-]\d\.\d\d)'
                        r'_r\+0\.00_s\+0\.00\.mod$')


def read_stars(path=HERE / 'stars.csv'):
    import pandas as pd
    return pd.read_csv(path, comment='#')


def av3d_params(name):
    m = AV3D_NAME.match(name)
    t = int(m['t'])
    teff = t if t > 1000 else 100 * t
    z = int(m['z']) / 10 * (-1 if m['s'] == 'm' else 1)
    return np.array([teff / 1000, int(m['g']) / 10, z])


def read_av3d(path):
    rows = Path(path).read_text().splitlines()
    data = np.array([[float(x) for x in r.split()] for r in rows[3:] if r.strip()])
    ltau, temp, pe, rho, vmic = data.T
    return dict(ltau=ltau, temp=temp, lpe=np.log(pe), lrho=np.log(rho), vmic=vmic)


def read_marcs(path):
    """lgTau5, T, ln Pe, ln rho of a MARCS .mod file (vmic is the model's constant)."""
    lines = Path(path).read_text(errors='ignore').split('\n')
    n = int(next(l for l in lines if 'Number of depth points' in l).split()[0])
    vt = float(next(l for l in lines if 'Microturbulence' in l).split()[0])
    i1 = next(i for i, l in enumerate(lines) if 'lgTauR' in l and 'Pe' in l)
    i2 = next(i for i, l in enumerate(lines) if 'lgTauR' in l and 'Density' in l)
    a = np.array([[float(x) for x in lines[i1 + 1 + k].split()[:7]] for k in range(n)])
    b = np.array([[float(x) for x in lines[i2 + 1 + k].split()[:4]] for k in range(n)])
    return dict(ltau=a[:, 2], temp=a[:, 4], lpe=np.log(a[:, 5]), lrho=np.log(b[:, 3]),
                vmic=np.full(n, vt))


def av3d_grid(folder):
    return {p.stem: (av3d_params(p.stem), p) for p in sorted(Path(folder).glob('*.atmos'))
            if AV3D_NAME.match(p.stem)}


def marcs_grid(folder, giant):
    """Standard-composition MARCS models of one geometry: {name: (params, path)}."""
    out = {}
    for p in Path(folder).glob('s*.mod' if giant else 'p*.mod'):
        m = MARCS_NAME.match(p.name)
        if not m:
            continue
        z = float(m['z'])
        alpha = min(0.4, max(0.0, -0.4 * z))
        std = (m['m'] == '1.0' and m['vt'] == '02') if giant else (m['m'] == '0.0' and m['vt'] == '01')
        if std and abs(float(m['a']) - round(alpha, 2)) < 0.01:
            out[p.stem] = (np.array([int(m['t']) / 1000, float(m['g']), z]), p)
    return out


def nearest(grid, x, k=K_NODES, exclude=()):
    names = [n for n in grid if n not in exclude]
    pts = np.array([grid[n][0] for n in names])
    d = np.sqrt((((pts - x) / SCALE) ** 2).sum(1))
    order = np.argsort(d)[:k]
    return [names[i] for i in order], d[order]


def interpolate(models, params, x, dist):
    """Weighted least-squares plane through ``models`` at ``params``, evaluated at ``x``,
    on the shared log tau range. Returns the model dict on NDEPTH points."""
    top = max(m['ltau'][0] for m in models)
    bottom = min(m['ltau'][-1] for m in models)
    ltau = np.linspace(top, bottom, NDEPTH)
    w = 1 / (dist ** 2 + 0.01)
    design = np.column_stack([np.ones(len(models)), (np.array(params) - x) / SCALE])
    rank = np.linalg.matrix_rank(design * np.sqrt(w)[:, None])
    if rank < design.shape[1]:
        # nodes in a plane (e.g. one metallicity): drop the degenerate directions
        keep = [0] + [j for j in range(1, 4) if np.ptp(design[:, j]) > 0]
        design = design[:, keep]
    out = {'ltau': ltau}
    for key in ('temp', 'lpe', 'lrho', 'vmic'):
        y = np.array([np.interp(ltau, m['ltau'], m[key]) for m in models])
        coef, *_ = np.linalg.lstsq(design * np.sqrt(w)[:, None], y * np.sqrt(w)[:, None], rcond=None)
        out[key] = coef[0]          # value at x (the plane's constant term)
    return out


def write_text(path, model, feh, name):
    rows = [name, f'{len(model["ltau"])} {feh}', '* depth T Pe rho vmic']
    rows += [f'{a:.6e} {t:.3f} {np.exp(p):.5e} {np.exp(r):.5e} {v:.4f}'
             for a, t, p, r, v in zip(model['ltau'], model['temp'], model['lpe'],
                                      model['lrho'], model['vmic'])]
    Path(path).write_text('\n'.join(rows) + '\n')


def build(av3d_dir, marcs_dir, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grids = {'av3d': av3d_grid(av3d_dir),
             'marcs_p': marcs_grid(marcs_dir, giant=False),
             'marcs_s': marcs_grid(marcs_dir, giant=True)}
    readers = {'av3d': read_av3d, 'marcs_p': read_marcs, 'marcs_s': read_marcs}
    log = ['# star kind nodes (scaled distance) | nearest node, offset dTeff dlogg dFeH']
    for s in read_stars().itertuples():
        x = np.array([s.teff / 1000, s.logg, s.feh])
        for kind in ('av3d', 'marcs_s' if s.logg < 3.5 else 'marcs_p'):
            names, dist = nearest(grids[kind], x)
            models = [readers[kind](grids[kind][n][1]) for n in names]
            if dist[0] < EXACT:
                names, dist, model = names[:1], dist[:1], models[0]
            else:
                model = interpolate(models, [grids[kind][n][0] for n in names], x, dist)
            label = 'av3d' if kind == 'av3d' else 'marcs'
            write_text(out_dir / f'{s.star}_{label}.atmos', model, s.feh, f'{s.star}_{label}')
            off = x - grids[kind][names[0]][0]
            log.append(f'{s.star} {label} ' + ' '.join(f'{n}({d:.2f})' for n, d in zip(names, dist))
                       + f' | {names[0]} {1000 * off[0]:+.0f} {off[1]:+.2f} {off[2]:+.2f}'
                       + f' | ltau {model["ltau"][0]:.2f}..{model["ltau"][-1]:.2f}')
    (out_dir / 'atmospheres.log').write_text('\n'.join(log) + '\n')
    print('\n'.join(log))


def validate(av3d_dir):
    """Leave-one-out error of the <3D> interpolant at every node: T, ln P_e and ln rho at
    log tau_500 = -3, -2, -1, 0 against the node's own average."""
    grid = av3d_grid(av3d_dir)
    probes = np.array([-3., -2., -1., 0.])
    rows = []
    for name, (x, path) in grid.items():
        names, dist = nearest(grid, x, exclude=(name,))
        if dist[0] > 1.5:
            continue
        model = interpolate([read_av3d(grid[n][1]) for n in names], [grid[n][0] for n in names], x, dist)
        truth = read_av3d(path)
        dt = np.interp(probes, model['ltau'], model['temp']) - np.interp(probes, truth['ltau'], truth['temp'])
        dr = np.interp(probes, model['ltau'], model['lrho']) - np.interp(probes, truth['ltau'], truth['lrho'])
        rows.append((name, dist[0], *dt, *dr))
    arr = np.array([r[1:] for r in rows])
    print(f'{len(rows)} nodes; nearest-neighbour distance median {np.median(arr[:, 0]):.2f}')
    for j, p in enumerate(probes):
        print(f'log tau {p:+.0f}: dT rms {np.sqrt(np.mean(arr[:, 1 + j] ** 2)):.0f} K, '
              f'p90 |dT| {np.percentile(np.abs(arr[:, 1 + j]), 90):.0f} K; '
              f'dln rho rms {np.sqrt(np.mean(arr[:, 5 + j] ** 2)):.3f}')
    return rows


if __name__ == '__main__':
    if sys.argv[1] == '--validate':
        validate(sys.argv[2])
    else:
        build(*sys.argv[1:4])
