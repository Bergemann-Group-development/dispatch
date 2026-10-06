"""Figures of the phase-6 report from the delta_maps.py / pilot_compare.py tables.

    python plot_phase6.py o O_p6.csv OUT.png          # Delta_NLTE of O I 6158 over the grid
    python plot_phase6.py pilot pilot_p6.csv OUT.png  # the pilot's paired terms per model
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_phase5 import style, INK, MUTED, COLOR, NAME

EW_MIN = 1.0     # mA; weaker O I 6158 is not measurable and its Delta is noise


def o_map(csv, out):
    d = pd.read_csv(csv)
    d = d[(d.line < 6200) & (d.ew_lte > EW_MIN)]
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, (title, sel) in zip(axs, (('dwarfs, log g >= 3.5', d.logg >= 3.5),
                                      ('giants, log g < 3.5', d.logg < 3.5))):
        g = d[sel]
        sc = ax.scatter(g.teff + 30 * (g.logg - g.logg.mean()), g.dA_nlte, c=g.feh, s=8,
                        cmap='viridis', vmin=-4, vmax=0.5)
        ax.set_title(title, fontsize=10, color=INK)
        ax.set_xlabel('Teff (K)', color=INK, fontsize=10)
        style(ax, 'Delta_NLTE O I 6158 (dex)' if ax is axs[0] else '')
    cb = fig.colorbar(sc, ax=axs, pad=0.01)
    cb.set_label('[Fe/H]', color=INK)
    fig.savefig(out, dpi=120, bbox_inches='tight')


def pilot(csv, out):
    """Median over each model's snapshots and lines of an element: D_nlte_3d, D_3d_lte_av,
    D_int_av and D_map, Mg I and Fe I side by side."""
    d = pd.read_csv(csv)
    terms = [t for t in ('D_nlte_3d', 'D_nlte_av', 'D_3d_lte_av', 'D_int_av', 'D_map') if t in d]
    models = list(dict.fromkeys(d.model))
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, el in zip(axs, ('Mg', 'Fe')):
        g = d[d.element == el].groupby('model', sort=False)[terms].median()
        x = np.arange(len(models))
        for i, t in enumerate(terms):
            y = [g.loc[m, t] if m in g.index else np.nan for m in models]
            ax.bar(x + (i - len(terms) / 2 + 0.5) * 0.16, y, 0.16, label=t)
        ax.set_xticks(x, models, rotation=40, fontsize=8, color=MUTED)
        ax.set_title(NAME[el], fontsize=10, color=COLOR[el])
        style(ax, 'dex' if ax is axs[0] else '')
    axs[0].legend(fontsize=7, frameon=False)
    fig.savefig(out, dpi=120, bbox_inches='tight')


if __name__ == '__main__':
    if len(sys.argv) != 4 or sys.argv[1] not in ('o', 'pilot'):
        sys.exit(__doc__)
    (o_map if sys.argv[1] == 'o' else pilot)(sys.argv[2], sys.argv[3])
