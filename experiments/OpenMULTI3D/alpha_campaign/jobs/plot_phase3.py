"""Figures and summary tables of phase 3, from the tables of delta_maps.py V and X.

    python plot_phase3.py vmic DELTA_V.csv OUT.png
    python plot_phase3.py lte3d DELTA_X.csv OUT.png
    python plot_phase3.py placement S4_LABELS.csv STAGGER3D.txt OUT.png

vmic: Delta_<3D> = A(<3D>) - A(MARCS) in LTE per model (median over the Fe I lines,
the Fe II lines and the Mg lines), with the <3D> average's own depth-dependent vmic and
with its MARCS twin's constant vmic, against Teff for giants (log g <= 3) and dwarfs.
lte3d: Delta_3D(LTE) = A(3D) - A(1D) per snapshot against the <3D> average (own vmic)
and the MARCS twin, per line group.  Both print the medians per [Fe/H].
placement: the HRS stars on their Gaia XP XGBoost labels (s4_labels.py) in the Kiel and
Teff-[M/H] planes with the 3D Stagger models, and the offsets of the phase-2 crude 1D LTE
fit from them.
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_batch2 import INK, MUTED, GRID, BLUE, stagger_labels

GROUPS = {'Fe': 'Fe I (8 lines)', 'FeII': 'Fe II (5 lines)', 'Mg': 'Mg I (3 lines)'}
FEH_COLOURS = {0.5: '#b52f2f', 0.0: '#ec8c8b', -0.5: '#e3b23c', -1.0: '#7fb069',
               -1.5: '#3a9e8f', -2.0: '#3987e5', -3.0: '#1c5cab', -4.0: '#0d366b'}


def colour(feh):
    keys = np.array(sorted(FEH_COLOURS))
    return FEH_COLOURS[keys[np.abs(keys - feh).argmin()]]


def medians(df, key, cols):
    return df.groupby([key, 'element'], as_index=False).agg(
        {**{c: 'first' for c in ('teff', 'logg', 'feh')}, **{c: 'median' for c in cols}})


def vmic(table, out):
    df = medians(pd.read_csv(table), 'av3d', ['dA_av', 'dA_avv'])
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.4), sharex=True, sharey=True,
                             constrained_layout=True)
    for i, (gname, sel) in enumerate((('giants, log g <= 3', df['logg'] <= 3.0),
                                      ('dwarfs, log g >= 3.5', df['logg'] >= 3.5))):
        for j, (el, title) in enumerate(GROUPS.items()):
            ax = axes[i, j]
            x = df[sel & (df['element'] == el)]
            c = [colour(f) for f in x['feh']]
            ax.scatter(x['teff'], x['dA_av'], s=14, facecolor='none', edgecolor=c, lw=0.8)
            ax.scatter(x['teff'], x['dA_avv'], s=10, color=c, lw=0)
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.grid(color=GRID, lw=0.5)
            ax.set_xlim(7150, 3850)
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(f'{gname}\nDelta_<3D> LTE [dex]')
            if i == 1:
                ax.set_xlabel('Teff [K]')
    for f, c in FEH_COLOURS.items():
        axes[0, 0].scatter([], [], color=c, s=10, label=f'{f:+.1f}')
    axes[0, 0].scatter([], [], facecolor='none', edgecolor=INK, s=14, label='<3D> own vmic')
    axes[0, 0].scatter([], [], color=INK, s=10, label='twin vmic')
    axes[0, 0].legend(fontsize=5.5, frameon=False, ncol=2, loc='upper left')
    fig.suptitle('Delta_<3D> = A(<3D>) - A(MARCS), LTE, median per line group', color=INK)
    fig.savefig(out)
    df['giant'] = df['logg'] <= 3.0
    print(df.groupby(['element', 'giant', 'feh'])[['dA_av', 'dA_avv']].median().round(3).to_string())
    cool = df['giant'] & (df['teff'] <= 4500) & (df['feh'] >= -0.5)
    print('cool metal-rich giants:\n', df[cool].groupby('element')[['dA_av', 'dA_avv']]
          .agg(['median', 'min', 'max']).round(3).to_string())
    print(f'-> {out}')


def lte3d(table, out):
    raw = pd.read_csv(table)
    cols = [c for c in ('dA_3d_av', 'dA_3d_avv', 'dA_3d_tw') if c in raw.columns]
    df = medians(raw, 'run', cols)
    refs = [(c, t) for c, t in (('dA_3d_av', 'vs <3D>'), ('dA_3d_tw', 'vs MARCS')) if c in df]
    fig, axes = plt.subplots(len(refs), 3, figsize=(7.4, 2.2 * len(refs) + 0.3), sharex=True,
                             sharey=True, constrained_layout=True, squeeze=False)
    for i, (c, rname) in enumerate(refs):
        for j, (el, title) in enumerate(GROUPS.items()):
            ax = axes[i, j]
            x = df[df['element'] == el]
            for giant, m in ((True, '^'), (False, 'o')):
                y = x[(x['logg'] <= 3.0) == giant]
                ax.scatter(y['teff'], y[c], s=14, marker=m, color=[colour(f) for f in y['feh']], lw=0)
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.grid(color=GRID, lw=0.5)
            ax.set_xlim(7150, 3850)
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(f'Delta_3D LTE {rname} [dex]')
            if i == len(refs) - 1:
                ax.set_xlabel('Teff [K]')
    for f, col in FEH_COLOURS.items():
        axes[0, 0].scatter([], [], color=col, s=10, label=f'{f:+.1f}')
    axes[0, 0].scatter([], [], marker='^', color=INK, s=10, label='log g <= 3')
    axes[0, 0].legend(fontsize=5.5, frameon=False, ncol=2, loc='upper left')
    fig.suptitle(f'Delta_3D(LTE) = A(3D) - A(1D), {df["run"].nunique()} snapshots, '
                 'median per line group', color=INK)
    fig.savefig(out)
    print(df.groupby(['element', 'feh'])[cols].agg(['median', 'count']).round(3).to_string())
    print(f'-> {out}')


def placement(table, stagger, out):
    d = pd.read_csv(table, dtype={'set': str})
    d['set'] = d['set'].str.zfill(4)
    lab = d[d['label_src'].fillna('') != '']
    s3 = stagger_labels([l.split()[0] for l in open(stagger) if l.strip()])
    near = np.array([((abs(s3[:, 0] - T) <= 250) & (abs(s3[:, 1] - g) <= 0.5)
                      & (abs(s3[:, 2] - m) <= 0.5)).any()
                     for T, g, m in zip(lab['teff'], lab['logg'], lab['mh'])])
    sets = {'0201': BLUE[4], '0301': '#e3b23c', '0401': '#7fb069', '0408': BLUE[2]}
    fig, axes = plt.subplots(1, 4, figsize=(9.6, 2.9), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.2, 1.2, 1, 1]))
    for ax, y, ylab, lim in zip(axes[:2], ('logg', 'mh'), ('log g', '[M/H]'), ((5.3, 0.7), (-4.3, 0.75))):
        ax.grid(color=GRID, lw=0.5)
        ys = s3[:, 1] if y == 'logg' else s3[:, 2]
        ax.scatter(s3[:, 0], ys, marker='s', s=30, facecolor='none', edgecolor=INK, lw=0.6,
                   zorder=3, label='3D Stagger model')
        for s, c in sets.items():
            x = lab[lab['set'] == s]
            ax.scatter(x['teff'], x[y], s=4, color=c, alpha=0.7, lw=0, label=f'{s} ({len(x)})')
        ax.set_xlim(7150, 3850)
        ax.set_ylim(*lim)
        ax.set_xlabel('Teff [K]')
        ax.set_ylabel(ylab)
    axes[0].legend(loc='upper left', frameon=False, fontsize=6, markerscale=2)
    c = lab[lab['teff_crude'].notna()]
    for ax, (q, cq, rng, unit) in zip(axes[2:], (('teff', 'teff_crude', (-1500, 1500), 'K'),
                                                 ('mh', 'feh_crude', (-3, 1.5), 'dex'))):
        for s, col in sets.items():
            dd = (c.loc[c['set'] == s, cq] - c.loc[c['set'] == s, q]).dropna()
            if len(dd):
                ax.hist(dd, bins=40, range=rng, histtype='step', color=col, label=s)
        ax.axvline(0, color=MUTED, lw=0.6)
        ax.set_xlabel(f'crude 1D LTE - XGBoost, {q} [{unit}]')
    axes[2].legend(frameon=False, fontsize=6)
    fig.suptitle(f'{len(lab)} HRS stars on Gaia XP XGBoost labels; {near.mean():.0%} within one '
                 'grid step of a 3D Stagger model', color=INK)
    fig.savefig(out)
    print(f'{len(lab)} stars, {near.mean():.3f} near a 3D model; by set: ' +
          ', '.join(f'{s}: {near[(lab["set"] == s).to_numpy()].mean():.2f}' for s in sets
                    if (lab['set'] == s).any()) + f' -> {out}')


if __name__ == '__main__':
    if sys.argv[1:2] == ['vmic']:
        vmic(*sys.argv[2:4])
    elif sys.argv[1:2] == ['placement']:
        placement(*sys.argv[2:5])
    elif sys.argv[1:2] == ['lte3d']:
        lte3d(*sys.argv[2:4])
    else:
        sys.exit(__doc__)
