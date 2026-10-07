#!/usr/bin/env python3
"""Figures of the atom-validation note (docs/notes/atom_validation) from the outputs of
``analyse.py`` and the extracted runs.

- ``fits_<el>.pdf``: observed line profiles against the NLTE syntheses of the original and
  the repaired atom at each atom's own abundance, with the EW window and the continuum
  level, for one star of each class; the lines are a rule-based sample of used, borderline,
  outlier, rejected and mis-matched lines.
- ``robust_window.pdf``, ``robust_instrument.pdf``: abundance change when the window edges,
  its position or the continuum level are perturbed, and the spectrum-to-spectrum offsets.
- ``fe_trends.pdf``: Fe I abundances against excitation energy, REW and wavelength.
- ``metrics.pdf``: per element and star, ionisation balance, slopes and [X/H] - GBS.
- ``clv_fe.pdf``, ``clv_others.pdf``: A(mu) - A(1) and dA/dmu for IAG, KPNO and Takeda.

    python plots.py LINES_DIR RUN_DIR SPECTRA_DIR ANA_DIR FIG_DIR [figure ...]
"""
import sys
from functools import lru_cache
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyse as A

COLOR = {'orig': '#eb6834', 'rep': '#2a78d6', 'mid': '#1baf7a', 'lte': '#8a8a8a', 'obs': '#1a1a1a',
         'window': '#e6e6e6'}
ELEMENT_COLOR = dict(zip(('Fe', 'Mg', 'Ca', 'Si', 'Al', 'Mn'),
                         ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300')))
STYLE = {'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
         'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6,
         'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'lines.linewidth': 1.0,
         'lines.markersize': 3, 'legend.frameon': False, 'pdf.fonttype': 42,
         'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02}
WIDTH = 6.3                                  # inches, the text width of docs/notes (16 cm)
CLASSES = (('Sun', 'Sun'), ('epsEri', 'dwarf'), ('HD140283', 'metal-poor dwarf'),
           ('Arcturus', 'giant'), ('HD122563', 'metal-poor giant'))
ORIGINAL = {el: pairs[0][0] for el, pairs in A.PAIRS.items()}
ELEMENTS = tuple(A.PAIRS)
ION = {1: 'I', 2: 'II'}


class Data:
    """The analysis tables and cached spectra the figures share."""

    def __init__(self, lines_dir, run_dir, spectra_dir, ana_dir):
        self.run_dir, self.spectra_dir, self.ana = Path(run_dir), Path(spectra_dir), Path(ana_dir)
        self.lines = pd.read_csv(Path(lines_dir) / 'lines.csv')
        self.res = A.flag_common(pd.read_csv(self.ana / 'line_abundances.csv'))
        self.runs = pd.read_csv(self.run_dir / 'runs_b.csv')
        self.stars = pd.read_csv(A.HERE / 'stars.csv', comment='#').set_index('star')

    @lru_cache(maxsize=None)
    def ews(self, star):
        return A.drop_discrepant(pd.read_csv(self.ana / f'ew_{star}.csv'))

    @lru_cache(maxsize=None)
    def obs(self, spectrum):
        return A.read_obs(self.spectra_dir / spectrum)

    @lru_cache(maxsize=None)
    def synth(self, star, atom):
        run = self.runs[(self.runs.star == star) & (self.runs.atom == atom) & (self.runs.atmos == 'av3d')].iloc[0]
        return A.load_synth(self.run_dir / 'extract' / 'b' / f'{run["name"]}.npz')


def status(r):
    """Role of one line result in the metrics."""
    if r.mismatch:
        return 'mismatch'
    if not r.clean:
        return 'rejected'
    if not r.common:
        return 'not common'
    if r.outlier:
        return 'outlier'
    borderline = r.blend > 0.2 or r.rew > -4.95 or r.w_obs < 8
    return 'borderline' if borderline else 'used'


def sample_lines(data, el):
    """A fixed rule picks the lines of the fit figure from the repaired atom's results in
    the class stars: the weakest and the strongest line used in most stars, the most often
    borderline line, a line rejected as an outlier, the most often rejected line, and up to
    two mis-matched lines."""
    stars = [s for s, _ in CLASSES]
    r = data.res[(data.res.element == el) & (data.res.atom == A.REPAIRED[el]) & (data.res.atmos == 'av3d')
                 & data.res.star.isin(stars)].copy()
    r['status'] = [status(x) for x in r.itertuples()]
    by = r.groupby(['ion', 'lam'])
    count = {k: by.status.apply(lambda s, k=k: (s == k).sum()) for k in
             ('used', 'borderline', 'outlier', 'rejected', 'mismatch')}
    w = by.w_obs.median()
    picks = []

    def add(index):
        for key in index:
            if key not in picks:
                picks.append(key)
                return

    used = count['used'][count['used'] == count['used'].max()].index
    add(w[used].sort_values().index)
    add(w[used].sort_values(ascending=False).index)
    for kind in ('borderline', 'outlier', 'rejected'):
        c = count[kind][count[kind] > 0]
        add(c.sort_values(ascending=False).index)
    for key in count['mismatch'][count['mismatch'] > 0].index[:2]:
        add([key])
    return picks, r.set_index(['star', 'ion', 'lam'])


def at_abundance(lam, spec_b, synth, da_fit):
    """NLTE ladder interpolated linearly at da_fit (da = 0 where the inversion failed)."""
    da, idx = A.ladder(synth, 'nlte')
    da_fit = 0.0 if not np.isfinite(da_fit) else np.clip(da_fit, da[0], da[-1])
    j = int(np.clip(np.searchsorted(da, da_fit) - 1, 0, len(da) - 2))
    t = (da_fit - da[j]) / (da[j + 1] - da[j])
    return (1 - t) * spec_b[idx[j]] + t * spec_b[idx[j + 1]]


def fit_panel(ax, data, star, el, ion, lam0, info):
    ews = data.ews(star)
    rows = ews[(ews.element == el) & (ews.ion == ion) & np.isclose(ews.lam, lam0) & (ews.atmos == 'av3d')]
    if rows.empty:
        ax.text(0.5, 0.5, 'not observed', transform=ax.transAxes, ha='center', va='center', fontsize=7)
        return
    best = rows[rows.atom == A.REPAIRED[el]].sort_values('snr').iloc[-1]
    lam_o, flux_o = data.obs(best.spectrum)
    lam_o = lam_o / (1 + best['shift'] / A.C_KMS)
    half = max(1.5 * (best.hi - best.lo), 0.6)
    sel = (lam_o > lam0 - half) & (lam_o < lam0 + half)
    ax.axvspan(best.lo, best.hi, color=COLOR['window'], lw=0)
    ax.axhline(best.norm, color=COLOR['lte'], lw=0.6, ls=':')
    ax.plot(lam_o[sel], flux_o[sel], color=COLOR['obs'], lw=0.6, drawstyle='steps-mid')
    fitted = {}
    for kind, atom in (('orig', ORIGINAL[el]), ('rep', A.REPAIRED[el])):
        syn = data.synth(star, atom)
        cut = (syn['lam'] > lam0 - half - 2) & (syn['lam'] < lam0 + half + 2)
        lam = syn['lam'][cut]
        spec_b = A.broaden(lam, syn['spec'][:, cut], best.fwhm)
        res_atom = data.res[(data.res.star == star) & (data.res.element == el) & (data.res.ion == ion)
                            & np.isclose(data.res.lam, lam0) & (data.res.atom == atom) & (data.res.atmos == 'av3d')]
        a = res_atom.a_nlte.iloc[0] if len(res_atom) else np.nan
        prof = at_abundance(lam, spec_b, syn, a - best.a_ref)
        keep = (lam > lam0 - half) & (lam < lam0 + half)
        ax.plot(lam[keep], best.norm * prof[keep], color=COLOR[kind], lw=0.9, ls='-' if kind == 'orig' else '--')
        fitted[kind] = a
    st = info.loc[(star, ion, lam0)] if (star, ion, lam0) in info.index else None
    label = status(st) if st is not None else ''
    da = fitted['rep'] - fitted['orig']
    text = f'{label}, W = {best.w_obs:.0f} mÅ'
    if st is not None and st.clean and np.isfinite(da):          # inversions of unclean lines are not meaningful
        text += f'\nΔA = {da:+.3f}'
    ax.text(0.03, 0.04, text, transform=ax.transAxes, fontsize=5.5, va='bottom', ha='left',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=0.8))
    ax.set_xlim(lam0 - half, lam0 + half)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.ticklabel_format(useOffset=False, axis='both')
    ax.tick_params(labelsize=6)


def fig_fits(data, el, fig_dir):
    picks, info = sample_lines(data, el)
    nrow, ncol = len(picks), len(CLASSES)
    fig, axes = plt.subplots(nrow, ncol, figsize=(WIDTH, 1.25 * nrow + 0.6), squeeze=False)
    for i, (ion, lam0) in enumerate(picks):
        for j, (star, cls) in enumerate(CLASSES):
            ax = axes[i, j]
            fit_panel(ax, data, star, el, ion, lam0, info)
            if i == 0:
                ax.set_title(f'{star} ({cls})')
            if j == 0:
                ax.set_ylabel(f'{el} {ION[ion]} {lam0:.2f}', fontsize=7)
    handles = [plt.Line2D([], [], color=COLOR['obs'], lw=0.6, label='observed'),
               plt.Line2D([], [], color=COLOR['orig'], label=ORIGINAL[el]),
               plt.Line2D([], [], color=COLOR['rep'], ls='--', label=A.REPAIRED[el]),
               plt.Line2D([], [], color=COLOR['lte'], ls=':', lw=0.6, label='continuum level'),
               plt.Rectangle((0, 0), 1, 1, color=COLOR['window'], label='EW window')]
    height = fig.get_figheight()
    fig.legend(handles=handles, loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.22 / height, 'air wavelength [Å]', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, 0.4 / height, 1, 1), h_pad=0.4, w_pad=0.3)
    fig.savefig(fig_dir / f'fits_{el}.pdf')
    plt.close(fig)


def fig_robust_window(data, fig_dir):
    t = pd.read_csv(data.ana / 'window_tests.csv')
    t = t[t.test != 'nominal']
    tests = (('edges', 'each edge moved outward [FWHM]'), ('shift', 'window shift [FWHM]'),
             ('continuum', 'continuum level change'))
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.3), sharey=True)
    for ax, (test, xlabel) in zip(axes, tests):
        g = t[t.test == test]
        steps = np.sort(g.step.unique())
        width = np.ptp(steps) / len(steps) / 8 if len(steps) > 1 else 0.01
        for k, el in enumerate(ELEMENTS):
            q = g[g.element == el].groupby('step').delta.quantile([0.16, 0.5, 0.84]).unstack()
            if q.empty:
                continue
            x = q.index.to_numpy() + (k - 2.5) * width
            ax.errorbar(x, q[0.5], yerr=[q[0.5] - q[0.16], q[0.84] - q[0.5]], fmt='o', ms=3, lw=0.8,
                        color=ELEMENT_COLOR[el], label=el)
        ax.axhline(0, color=COLOR['lte'], lw=0.5)
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel('ΔA(NLTE) [dex], median and 16–84%')
    axes[0].legend(ncol=2, loc='upper left')
    fig.tight_layout(w_pad=0.5)
    fig.savefig(fig_dir / 'robust_window.pdf')
    plt.close(fig)


def fig_robust_instrument(data, fig_dir):
    s = pd.read_csv(data.ana / 'spectrum_abundances.csv')
    s = s[(s.element == 'Fe') & (s.ion == 1)]
    order = sorted(s.star.unique(), key=lambda x: -data.stars.loc[x, 'teff'])
    instruments = ('HARPS', 'NARVAL', 'UVES')
    inst_color = dict(zip(instruments, ('#2a78d6', '#eb6834', '#1baf7a')))
    fig, ax = plt.subplots(figsize=(WIDTH, 2.2))
    for i, star in enumerate(order):
        g = s[s.star == star]
        files = sorted(g.spectrum.unique())
        for k, f in enumerate(files):
            q = g[g.spectrum == f].offset.quantile([0.16, 0.5, 0.84])
            inst = f.split('/')[0]
            x = i + (k - (len(files) - 1) / 2) * 0.12
            ax.errorbar(x, q[0.5], yerr=[[q[0.5] - q[0.16]], [q[0.84] - q[0.5]]], fmt='o', ms=3, lw=0.8,
                        color=inst_color.get(inst, COLOR['lte']))
    for inst in instruments:
        ax.plot([], [], 'o', color=inst_color[inst], label=inst)
    ax.axhline(0, color=COLOR['lte'], lw=0.5)
    ax.set_xticks(range(len(order)), order, rotation=40, ha='right')
    ax.set_ylabel('A(Fe I) − mean over spectra [dex]')
    ax.legend(ncol=3, loc='upper right')
    fig.tight_layout()
    fig.savefig(fig_dir / 'robust_instrument.pdf')
    plt.close(fig)


def _fit_line(ax, x, y, color, ls):
    if len(x) >= 4 and np.ptp(x) > 0:
        c = np.polyfit(x, y, 1)
        xx = np.array([x.min(), x.max()])
        ax.plot(xx, np.polyval(c, xx), color=color, ls=ls, lw=1)
        return c[0]
    return np.nan


def fig_fe_trends(data, fig_dir):
    use = pd.read_csv(data.ana / 'line_abundances_used.csv')
    use = use[(use.element == 'Fe') & (use.atmos == 'av3d')]
    cols = (('evi', 'χ [eV]', 1.0), ('rew', 'log(W/λ)', 1.0), ('lam', 'λ [1000 Å]', 1e-3))
    fig, axes = plt.subplots(len(CLASSES), 3, figsize=(WIDTH, 1.5 * len(CLASSES) + 0.3), sharex='col')
    for i, (star, cls) in enumerate(CLASSES):
        g = use[use.star == star]
        for j, (col, label, scale) in enumerate(cols):
            ax = axes[i, j]
            text = []
            for kind, atom, ls in (('orig', 'fe607a_pd', '-'), ('rep', 'fe607a_fix', '--')):
                a = g[g.atom == atom]
                one, two = a[a.ion == 1], a[a.ion == 2]
                ax.plot(one[col] * scale, one.a_nlte, '.', ms=2, color=COLOR[kind], alpha=0.5)
                ax.plot(two[col] * scale, two.a_nlte, 's', ms=2.5, mfc='none', mew=0.6, color=COLOR[kind])
                s = _fit_line(ax, one[col].to_numpy() * scale, one.a_nlte.to_numpy(), COLOR[kind], ls)
                text.append(f'{s:+.3f}')
            ax.text(0.98, 0.04, ' / '.join(text), transform=ax.transAxes, ha='right', fontsize=6)
            if j == 0:
                ax.set_ylabel(f'{star}\nA(Fe)', fontsize=7)
            if i == len(CLASSES) - 1:
                ax.set_xlabel(label)
            lo, hi = np.nanpercentile(g.a_nlte, [1, 99]) if len(g) else (0, 1)
            ax.set_ylim(lo - 0.1, hi + 0.1)
    handles = [plt.Line2D([], [], color=COLOR['orig'], marker='.', label='fe607a_pd, Fe I'),
               plt.Line2D([], [], color=COLOR['rep'], marker='.', ls='--', label='fe607a_fix, Fe I'),
               plt.Line2D([], [], color=COLOR['lte'], marker='s', mfc='none', ls='none', label='Fe II')]
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.03, 1, 1), h_pad=0.3)
    fig.savefig(fig_dir / 'fe_trends.pdf')
    plt.close(fig)


def fig_metrics(data, fig_dir):
    met = pd.read_csv(data.ana / 'metrics.csv')
    met = met[met.atmos == 'av3d']
    s2s = pd.read_csv(data.ana / 'star_to_star.csv')
    s2s = s2s[(s2s.atmos == 'av3d') & (s2s.star != 'Sun')]
    rows = (('ion_bal', 'ion_bal_err', 'A(II) − A(I)'), ('chi_slope_I', 'chi_slope_err_I', 'dA/dχ [dex/eV]'),
            ('rew_slope_I', 'rew_slope_err_I', 'dA/dlog(W/λ)'), ('lam_slope_I', 'lam_slope_err_I', 'dA/dλ [dex/kÅ]'),
            ('d_gbs', 'e_gbs', '[X/H] − GBS'))
    fig, axes = plt.subplots(len(rows), len(ELEMENTS), figsize=(WIDTH, 7.5), sharex=True)
    for j, el in enumerate(ELEMENTS):
        for i, (col, err, label) in enumerate(rows):
            ax = axes[i, j]
            src = s2s if col == 'd_gbs' else met
            drawn = False
            for atom, kind, dx in ((ORIGINAL[el], 'orig', -40), (A.REPAIRED[el], 'rep', 40)):
                g = src[(src.element == el) & (src.atom == atom)]
                if col not in g or g[col].isna().all():
                    continue
                g = g.merge(data.stars[['teff']], left_on='star', right_index=True, suffixes=('_', ''))
                e = g[err] if err in g and col != 'd_gbs' else None
                ax.errorbar(g.teff + dx, g[col], yerr=e, fmt='o', ms=2.5, lw=0.6, color=COLOR[kind])
                drawn = True
            ax.axhline(0, color=COLOR['lte'], lw=0.5)
            if not drawn:
                ax.text(0.5, 0.5, 'n/a', transform=ax.transAxes, ha='center', va='center', fontsize=7)
                ax.set_yticks([])
            if i == 0:
                ax.set_title(el)
            if j == 0:
                ax.set_ylabel(label, fontsize=7)
            if i == len(rows) - 1:
                ax.set_xlabel('Teff [K]')
    axes[-1, 0].invert_xaxis()
    handles = [plt.Line2D([], [], color=COLOR['orig'], marker='o', ls='none', label='original'),
               plt.Line2D([], [], color=COLOR['rep'], marker='o', ls='none', label='repaired')]
    fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.02, 1, 1), h_pad=0.3, w_pad=0.3)
    fig.savefig(fig_dir / 'metrics.pdf')
    plt.close(fig)


MU_BINS = np.array([0.15, 0.25, 0.32, 0.37, 0.45, 0.55, 0.65, 0.75, 0.85, 0.925, 0.96, 0.975, 0.985, 0.995, 1.01])


def _relative(long):
    """A(mu) - A(mu = 1) per line, atom and kind (Takeda: mu binned to MU_BINS)."""
    long = long.copy()
    long['mu_bin'] = pd.cut(long.mu, MU_BINS, labels=False)
    keys = A.CLV_KEYS + ['kind']
    long = long.groupby(keys + ['mu_bin'], as_index=False).agg(mu=('mu', 'mean'), a=('a', 'mean'))
    centre = long[long.mu >= 0.99].groupby(keys).a.mean()
    long['da'] = long.a - centre.reindex(pd.MultiIndex.from_frame(long[keys])).to_numpy()
    return long.dropna(subset=['da'])


def fig_clv_fe(data, fig_dir):
    long = pd.read_csv(data.ana / 'clv_common.csv')
    long = _relative(long[(long.element == 'Fe') & (long.atmos == 'av3d') & (long.ion == 1)])
    slopes = pd.read_csv(data.ana / 'clv_slopes.csv')
    slopes = slopes[(slopes.element == 'Fe') & (slopes.atmos == 'av3d') & (slopes.ion == 1)]
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 4.2), sharey='row')
    curves = (('lte', 'fe607a_pd', 'lte', 'LTE'), ('nlte', 'fe607a_pd', 'orig', 'fe607a_pd'),
              ('nlte', 'fe607a_fix', 'rep', 'fe607a_fix'))
    for j, src in enumerate(('IAG', 'KPNO', 'Takeda')):
        ax = axes[0, j]
        for kind, atom, color, label in curves:
            g = long[(long.source == src) & (long.kind == kind) & (long.atom == atom)]
            q = g.groupby('mu_bin').agg(mu=('mu', 'mean'), lo=('da', lambda x: x.quantile(0.16)),
                                        med=('da', 'median'), hi=('da', lambda x: x.quantile(0.84)))
            ax.fill_between(q.mu, q.lo, q.hi, color=COLOR[color], alpha=0.15, lw=0)
            ax.plot(q.mu, q.med, 'o-', ms=2.5, color=COLOR[color], label=label)
        n = long[(long.source == src)][['ion', 'lam']].drop_duplicates().shape[0]
        ax.set_title(f'{src} ({n} Fe I lines)')
        ax.axhline(0, color=COLOR['lte'], lw=0.5)
        ax.set_xlabel('μ')
        ax.invert_xaxis()
        ax2 = axes[1, j]
        bins = np.linspace(-1.2, 0.8, 41)
        for kind, atom, color, label in curves:
            g = slopes[(slopes.source == src) & (slopes.kind == kind) & (slopes.atom == atom)]
            ax2.hist(g.slope, bins=bins, histtype='step', color=COLOR[color], lw=1,
                     label=f'{label}: {g.slope.mean():+.3f}')
        ax2.axvline(0, color=COLOR['lte'], lw=0.5)
        ax2.set_xlabel('dA/dμ [dex]')
        ax2.legend(loc='upper left', fontsize=6)
    axes[0, 0].set_ylabel('A(μ) − A(1) [dex]')
    axes[1, 0].set_ylabel('lines')
    axes[0, 0].legend(loc='lower left')
    fig.tight_layout(h_pad=0.6)
    fig.savefig(fig_dir / 'clv_fe.pdf')
    plt.close(fig)


def fig_clv_others(data, fig_dir):
    long = pd.read_csv(data.ana / 'clv_common.csv')
    long = _relative(long[(long.element != 'Fe') & (long.atmos == 'av3d') & (long.kind == 'nlte')])
    elements = [el for el in ELEMENTS if el != 'Fe']
    sources = ('IAG', 'KPNO', 'Takeda')
    fig, axes = plt.subplots(len(sources), len(elements), figsize=(WIDTH, 5.0), sharex=True, sharey=True)
    for j, el in enumerate(elements):
        pairs = [(ORIGINAL[el], 'orig', '-'), (A.REPAIRED[el], 'rep', '--')]
        if el == 'Mg':
            pairs.insert(1, ('mg86e_pd', 'mid', ':'))
        for i, src in enumerate(sources):
            ax = axes[i, j]
            g = long[(long.element == el) & (long.source == src)]
            for atom, kind, ls in pairs:
                for _, line in g[g.atom == atom].groupby(['ion', 'lam']):
                    line = line.sort_values('mu')
                    ax.plot(line.mu, line.da, ls=ls, marker='o' if src == 'KPNO' else None, ms=2,
                            color=COLOR[kind], lw=0.7, alpha=0.8)
            if g.empty:
                ax.text(0.5, 0.5, 'no common lines', transform=ax.transAxes, ha='center', va='center', fontsize=6)
            ax.axhline(0, color=COLOR['lte'], lw=0.5)
            if i == 0:
                ax.set_title(el)
            if j == 0:
                ax.set_ylabel(f'{src}\nA(μ) − A(1)', fontsize=7)
            if i == len(sources) - 1:
                ax.set_xlabel('μ')
            ax.set_xticks([1.0, 0.6, 0.2])
    axes[0, 0].invert_xaxis()
    handles = [plt.Line2D([], [], color=COLOR['orig'], label='original'),
               plt.Line2D([], [], color=COLOR['mid'], ls=':', label='mg86e_pd'),
               plt.Line2D([], [], color=COLOR['rep'], ls='--', label='repaired')]
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.03, 1, 1), h_pad=0.3, w_pad=0.3)
    fig.savefig(fig_dir / 'clv_others.pdf')
    plt.close(fig)


FIGURES = {'fits': lambda d, f: [fig_fits(d, el, f) for el in ELEMENTS],
           'robust_window': fig_robust_window, 'robust_instrument': fig_robust_instrument,
           'fe_trends': fig_fe_trends, 'metrics': fig_metrics, 'clv_fe': fig_clv_fe,
           'clv_others': fig_clv_others}


def main(lines_dir, run_dir, spectra_dir, ana_dir, fig_dir, *names):
    plt.rcParams.update(STYLE)
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    data = Data(lines_dir, run_dir, spectra_dir, ana_dir)
    for name in names or FIGURES:
        FIGURES[name](data, fig_dir)
        print(f'{name}: done')


if __name__ == '__main__':
    main(*sys.argv[1:])
