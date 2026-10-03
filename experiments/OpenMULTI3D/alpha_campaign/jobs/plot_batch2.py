"""First-look figures of batch 2, from the tables the cluster jobs write.

    python plot_batch2.py placement LABELS.csv STAGGER3D.txt OUT.png
    python plot_batch2.py dnlte DELTA_A.csv OUT.png
    python plot_batch2.py d3d DELTA_B.csv OUT.png

placement: the 4MOST HRS stars (one point per star, its highest-S/N spectrum) in the
Kiel and Teff-[Fe/H] planes, on the crude 1D LTE labels of place_stars.py (log g held
at the Gaia value where there is one), with the 3D Stagger models (STAGGER3D.txt: one
folder name per line, e.g. t45g25m05_new) as open squares.
dnlte / d3d: abundance-like corrections per line across Teff and log g, in [Fe/H]
panels.
"""
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

INK, MUTED, GRID = '#262624', '#73726c', '#e2e1dc'
BLUE = ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
DIVERGING = LinearSegmentedColormap.from_list(
    'blue_red', ['#184f95', '#6da7ec', '#f0efec', '#ec8c8b', '#b52f2f'])
SET_NAMES = {201: '0201', 301: '0301', 401: '0401', 408: '0408'}
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.titlesize': 8.5,
                     'axes.titlecolor': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})


def stagger_labels(names):
    out = []
    for n in names:
        m = re.match(r't(\d+)g(\d\d)([mp])(\d\d)', n)
        teff = int(m[1]) * (100 if len(m[1]) == 2 else 1)
        out.append((teff, int(m[2]) / 10, int(m[4]) / 10 * (-1 if m[3] == 'm' else 1)))
    return np.unique(np.array(out), axis=0)


def placement(labels, stagger, out):
    d = pd.read_csv(labels, dtype={'specid': str}).drop_duplicates('specid')
    d['star'] = d['OBJECT'].fillna(d['specid'])
    d = d.sort_values('snr').drop_duplicates('star', keep='last')
    held = d['teff_c'].notna()
    d['T'] = np.where(held, d['teff_c'], d['teff'])
    d['g'] = np.where(held, d['logg_c'], d['logg'])
    d['m'] = np.where(held, d['feh_c'], d['feh'])
    s3 = stagger_labels([l.split()[0] for l in open(stagger) if l.strip()])
    near = np.array([((abs(s3[:, 0] - T) <= 250) & (abs(s3[:, 1] - g) <= 0.5)
                      & (abs(s3[:, 2] - m) <= 0.5)).any() for T, g, m in zip(d['T'], d['g'], d['m'])])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), constrained_layout=True)
    for ax, y, ylab, lim in zip(axes, ('g', 'm'), ('log g', '[Fe/H]'), ((5.3, 0.7), (-4.3, 0.75))):
        ax.grid(color=GRID, lw=0.5)
        ys = s3[:, 1] if y == 'g' else s3[:, 2]
        ax.scatter(s3[:, 0], ys, marker='s', s=40, facecolor='none', edgecolor=INK, lw=0.7,
                   zorder=3, label='3D Stagger model')
        n = d.groupby([d['T'].round(-1), d[y].round(2)]).size().reset_index(name='n')
        ax.scatter(n['T'], n[y], s=4 + 3 * n['n'], color=BLUE[2], alpha=0.75, lw=0, zorder=2,
                   label='4MOST stars (area = count)')
        for k in (1, 10, 50):
            ax.scatter([], [], s=4 + 3 * k, color=BLUE[2], alpha=0.75, lw=0, label=f'{k}')
        ax.set_xlim(7150, 3850)
        ax.set_ylim(*lim)
        ax.set_xlabel('Teff [K]')
        ax.set_ylabel(ylab)
    axes[0].legend(loc='upper left', frameon=False, fontsize=6.5, ncol=2)
    fig.suptitle(f'{len(d)} 4MOST HRS spectra sources on crude 1D LTE labels (not science labels); '
                 f'{near.mean():.0%} have a 3D model within one grid step', color=INK, fontsize=8)
    fig.savefig(out)
    print(f'{len(d)} stars, {held.sum()} with log g held at Gaia, {near.mean():.3f} near a 3D model; '
          f'by set: ' + ', '.join(f'{k}: {near[(d["set"] == k).to_numpy()].mean():.2f}'
                                  for k in sorted(d['set'].unique())) + f' -> {out}')


def panels(df, cols, out, title, label):
    """One row per [Fe/H] slice, one column per (title, column, row filter)."""
    fehs = [f for f in (0.0, -1.0, -2.0, -3.0) if (abs(df['feh'] - f) < 0.01).any()]
    fig, axes = plt.subplots(len(fehs), len(cols), figsize=(1.9 * len(cols) + 0.9, 1.75 * len(fehs)),
                             sharex=True, sharey=True, constrained_layout=True, squeeze=False)
    vals = np.concatenate([df.loc[sel, c].to_numpy() for _, c, sel in cols])
    lim = np.nanpercentile(np.abs(vals), 98) or 0.05
    norm = TwoSlopeNorm(0, -lim, lim)
    for i, f in enumerate(fehs):
        for j, (name, c, sel) in enumerate(cols):
            ax = axes[i, j]
            x = df[sel & (abs(df['feh'] - f) < 0.01)]
            sc = ax.scatter(x['teff'], x['logg'], c=x[c], cmap=DIVERGING, norm=norm, s=22,
                            marker='s', lw=0)
            ax.set_xlim(7150, 3850)
            ax.set_ylim(5.3, 0.7)
            ax.grid(color=GRID, lw=0.5)
            if i == 0:
                ax.set_title(name)
            if j == 0:
                ax.set_ylabel(f'[Fe/H] = {f:+.0f}\nlog g')
            if i == len(fehs) - 1:
                ax.set_xlabel('Teff [K]')
    fig.colorbar(sc, ax=axes, label=label, shrink=0.6)
    fig.suptitle(title, color=INK, fontsize=8.5)
    fig.savefig(out)
    print(f'-> {out}')


def line_median(df, col, element):
    """Median of col over the lines of element per model, as one pseudo-line."""
    key = [c for c in ('model', 'av3d') if c in df.columns][0]
    out = df[df['element'] == element].groupby(key, as_index=False).agg(
        {**{c: 'first' for c in ('teff', 'logg', 'feh')}, col: 'median'})
    return out.assign(element=f'{element} median', line=0.0)


def dnlte(table, out):
    df = pd.read_csv(table)
    df = pd.concat([df] + [line_median(df, 'dA_nlte', e) for e in ('Fe', 'FeII')
                           if (df['element'] == e).any()])
    cols = [(f'Mg I {l:.0f}', 'dA_nlte', (df['element'] == 'Mg') & (abs(df['line'] - l) < 0.1))
            for l in (5172.68, 5528.41, 5711.08)]
    cols += [('Fe I, median of 8', 'dA_nlte', df['element'] == 'Fe median')]
    if (df['element'] == 'FeII median').any():
        cols += [('Fe II, median of 5', 'dA_nlte', df['element'] == 'FeII median')]
    panels(df, cols, out, 'Delta_NLTE in 1D (MARCS), EW-based, A(NLTE) - A(LTE)',
           'abundance correction [dex]')


def d3d(table, out):
    df = pd.read_csv(table)
    extra = []
    for e in ('Fe', 'FeII'):
        if (df['element'] == e).any():
            lte, nlte = line_median(df, 'dA_3d_lte', e), line_median(df, 'dA_3d_nlte', e)
            extra.append(lte.merge(nlte[['av3d', 'dA_3d_nlte']], on='av3d'))
    df = pd.concat([df] + extra)
    mg = (df['element'] == 'Mg') & (abs(df['line'] - 5528.41) < 0.1)
    cols = [('Mg I 5528, LTE', 'dA_3d_lte', mg), ('Mg I 5528, NLTE', 'dA_3d_nlte', mg),
            ('Fe I median, LTE', 'dA_3d_lte', df['element'] == 'Fe median'),
            ('Fe I median, NLTE', 'dA_3d_nlte', df['element'] == 'Fe median')]
    if (df['element'] == 'FeII median').any():
        cols += [('Fe II median, NLTE', 'dA_3d_nlte', df['element'] == 'FeII median')]
    panels(df, cols, out, 'Delta_<3D>: A(<3D>) - A(MARCS) at the same EW and composition',
           'abundance correction [dex]')


if __name__ == '__main__':
    what = sys.argv[1:2]
    if what == ['placement']:
        placement(*sys.argv[2:5])
    elif what == ['dnlte']:
        dnlte(*sys.argv[2:4])
    elif what == ['d3d']:
        d3d(*sys.argv[2:4])
    else:
        sys.exit(__doc__)
