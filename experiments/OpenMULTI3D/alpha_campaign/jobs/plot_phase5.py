"""Figures of the phase-5 Delta summary from the delta_maps.py tables.

    python plot_phase5.py DELTA_DIR OLD_A2.csv OUT_DIR

DELTA_DIR holds A2_p5.csv, N_p5.csv and W_p5.csv (runs/delta/); OLD_A2.csv is the A2
table before the batch-2 redo (phase 4).  Writes delta_nlte_1d.png, delta_av3d.png and
delta_3d_vmic.png.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED = '#0b0b0b', '#52514e'
COLOR = {'Fe': BLUE, 'Mg': ORANGE, 'FeII': AQUA}
NAME = {'Fe': 'Fe I', 'Mg': 'Mg I', 'FeII': 'Fe II'}


def style(ax, ylabel):
    ax.axhline(0, color='#c3c2b7', lw=1, zorder=0)
    ax.grid(axis='y', color='#ecebe8', lw=0.8)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.set_ylabel(ylabel, color=INK, fontsize=10)


def nlte_1d(a2, old, out):
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for el in ('Fe', 'Mg'):
        for df, ls, tag in ((a2, '-', 'phase-4 setup'), (old, '--', 'batch 2')):
            d = df[(df.element == el) & ~df.get('emission', pd.Series(False, index=df.index))].groupby('feh')['dA_nlte'].median()
            ax.plot(d.index, d.values, ls, lw=2, color=COLOR[el])
            ax.annotate(f'{NAME[el]}, {tag}', (d.index[0], d.values[0]), xytext=(-4, 0),
                        textcoords='offset points', ha='right', va='center', fontsize=8, color=INK)
    ax.set_xlim(-5.4, 0.6)
    ax.set_xlabel('[Fe/H] of the MARCS model', color=INK, fontsize=10)
    style(ax, 'median Delta_NLTE (dex)')
    ax.set_title('1D Delta_NLTE, stage-2 SEDs: redo at [Fe/H] <= -1 (solid) vs batch 2 (dashed)',
                 fontsize=9, color=INK, loc='left')
    fig.tight_layout()
    fig.savefig(out / 'delta_nlte_1d.png', dpi=150)


def av3d(n, out):
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for el in ('Fe', 'Mg', 'FeII'):
        for col, ls, tag in (('dA_3d_lte_new', '-', 'LTE'), ('dA_3d_nlte_new', '--', 'NLTE')):
            if el == 'FeII' and tag == 'NLTE':
                continue
            d = n[n.element == el].groupby('feh')[col].median().dropna()
            ax.plot(d.index, d.values, ls, lw=2, color=COLOR[el], marker='o', ms=4)
            dy = {('Fe', 'NLTE'): 5, ('Mg', 'NLTE'): -5}.get((el, tag), 0)
            ax.annotate(f'{NAME[el]} {tag}', (d.index[0], d.values[0]), xytext=(-6, dy),
                        textcoords='offset points', ha='right', va='center', fontsize=8, color=INK)
    ax.set_xlim(-5.0, 0.7)
    ax.set_xlabel('[Fe/H] of the <3D> average', color=INK, fontsize=10)
    style(ax, 'median Delta_<3D> (dex)')
    ax.set_title('Own <3D> averages against their MARCS twins, A(<3D>) - A(MARCS)',
                 fontsize=9, color=INK, loc='left')
    fig.tight_layout()
    fig.savefig(out / 'delta_av3d.png', dpi=150)


def vmic(w, out):
    fe = w[w.element == 'Fe'].groupby('av3d')[['dA_3d_marcs_conv', 'dA_3d_marcs_vfit', 'vfit_marcs', 'teff', 'logg']].median().dropna()
    fe = fe.sort_values(['logg', 'teff'])
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    y = range(len(fe))
    ax.hlines(list(y), fe.dA_3d_marcs_conv, fe.dA_3d_marcs_vfit, color='#c3c2b7', lw=1.5)
    ax.plot(fe.dA_3d_marcs_conv, y, 'o', ms=6, color=ORANGE, label='vmic 1 / 2 km/s (MARCS convention)')
    ax.plot(fe.dA_3d_marcs_vfit, y, 'o', ms=6, color=BLUE, label='fitted vmic')
    ax.set_yticks(list(y))
    ax.set_yticklabels([f'{m}  ({v:.2f})' for m, v in zip(fe.index, fe.vfit_marcs)], fontsize=8)
    ax.axvline(0, color='#c3c2b7', lw=1, zorder=0)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel('median Delta_3D of the Fe I lines, A(3D) - A(MARCS) (dex)', color=INK, fontsize=9)
    fig.set_size_inches(6.4, 4.3)
    ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.4, -0.13), ncol=2)
    ax.set_title('3D LTE against MARCS; fitted vmic (km/s) in brackets', fontsize=9, color=INK, loc='left')
    fig.tight_layout()
    fig.savefig(out / 'delta_3d_vmic.png', dpi=150)


def main(delta_dir, old_a2, out):
    delta_dir, out = Path(delta_dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    nlte_1d(pd.read_csv(delta_dir / 'A2_p5.csv'), pd.read_csv(old_a2), out)
    av3d(pd.read_csv(delta_dir / 'N_p5.csv'), out)
    vmic(pd.read_csv(delta_dir / 'W_p5.csv'), out)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
