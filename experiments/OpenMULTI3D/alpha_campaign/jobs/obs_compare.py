"""The pilot's spectra against 4MOST HRS stars with labels near each pilot model.

    python obs_compare.py extract OUT.pkl            # cluster: line windows of the pilot runs
    python obs_compare.py compare WINDOWS.pkl OBS_DIR OUT.csv [TEFF LOGG MH]   # anywhere

extract reads, per pilot model, the normalised flux within +-WIN A of every LINES_SED line
of these spectra (pilot_compare.py has the same runs):
    3d_nlte   the snapshot's stage-2 SED on its own 3D departures (_3dmg4fe, _3dmgfe,
              _3dmg4, _3dmg, the first that exists; fe=True if it has Fe departures)
    3d_lte    the plain 3D LTE run (s<snap>_lte1, else _lte), since the SEDs' LTE set
              lacks the NLTE species' lines (reports/2026-10-05_phase5.md, section 5)
    map_nlte  the SED on the own <3D> departures mapped onto the snapshot (_map4mg[fe])
    av_nlte, av_lte   the own <3D> average's SED (b4/sed/av3d_<model>)
    tw_nlte, tw_lte   the MARCS twin's SED (b6/sed/<twin> at [Fe/H] <= -1, else p1sed/<model>/marcs)
and the curve-of-growth slopes dlog(EW)/dA of the nearest batch-1a run (delta_maps.slopes).

compare takes the stars whose S4 labels (s4_labels.csv) lie within TEFF/LOGG/MH (150 K,
0.3, 0.3) of a model, the highest-S/N spectrum of each (place_stars labels.csv, S/N >= 30),
from OBS_DIR/<set>/<stack>/.  Every model window is convolved to R = 20000 (Gaussian in
ln lambda) and put on the observed wavelengths.  Per star, line and spectrum the observed
flux is divided by a straight line fitted to observed/model outside the line core
(|lambda - lambda0| > hw, within WIN), as place_stars.py does for the whole arm; then

    dA = (log EW_obs - log EW_model) / slope

over the core (+-max(hw, HW_MIN)), the abundance change that would make the model match
the star.  3D spectra are averaged over the model's snapshots before the EW (3D NLTE
Fe I over the snapshots with Fe departures only).  The two
CSV files (labels.csv, s4_labels.csv) are read from the folder of WINDOWS.pkl.
Fe II rows keep their LTE value in the NLTE columns (delta_maps.flag_rows).
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WIN = 3.0           # A, half width of the stored window and of the continuum fit
HW_MIN = 0.4        # A, smallest core half width at R = 20000 (FWHM 0.26-0.33 A)
R_INST = 20000
KINDS = ('3d_nlte', '3d_lte', 'map_nlte', 'av_nlte', 'av_lte', 'tw_nlte', 'tw_lte')


def windows(run, lines, LTE=False):
    out = {}
    for _, lam0, _ in lines:
        s = run.crop(xmin=lam0 - WIN, xmax=lam0 + WIN, LTE=LTE)
        out[lam0] = (np.asarray(s.xdata.value, float), np.asarray(s.ydata.value, float))
    return out


def extract(out):
    import tumult
    from delta_maps import RUNS, LINES_SED, slopes, nearest_b1
    from make_nml import PILOT, run_name

    done = lambda f: (f / 'save').is_dir() and any((f / 'save').glob('flux_*.bin'))
    names = sorted(p.name for p in (RUNS / 'b1/hrs').iterdir())
    res = {}
    for folder, marcs, av, feh, alpha in PILOT:
        model = Path(folder).name.replace('_new', '')
        top = RUNS / 'p1sed' / model
        teff = 5777 if model.startswith('t5777') else int(model[1:3]) * 100
        logg = int(model.split('g')[1][:2]) / 10
        m = dict(teff=teff, logg=logg, feh=feh, alpha=alpha, spectra={}, snaps={},
                 slope=slopes(f'b1/hrs/{nearest_b1(teff, logg, feh, names)}', LINES_SED))
        tw = RUNS / f'b6/sed/{run_name(Path(marcs))}' if feh <= -1 else top / 'marcs'
        for kind, f in (('av', RUNS / f'b4/sed/av3d_{model}'), ('tw', tw)):
            if done(f):
                r = tumult.read(str(f))
                m['spectra'][f'{kind}_nlte'] = windows(r, LINES_SED)
                m['spectra'][f'{kind}_lte'] = windows(r, LINES_SED, LTE=True)
        for snap in sorted({p.name.split('_')[0] for p in top.glob('s*_3dmg*')}):
            found = [f for f in (top / f'{snap}_{k}' for k in ('3dmg4fe', '3dmgfe', '3dmg4', '3dmg'))
                     if done(f)]
            lte = [f for f in (RUNS / f'p1/{model}/{snap}_{k}' for k in ('lte1', 'lte')) if done(f)]
            if not found or not lte:
                continue
            fe = found[0].name.endswith('fe')
            sp = {'3d_nlte': windows(tumult.read(str(found[0])), LINES_SED),
                  '3d_lte': windows(tumult.read(str(lte[0])), LINES_SED)}
            mapped = top / f'{snap}_map4mg{"fe" if fe else ""}'
            if done(mapped):
                sp['map_nlte'] = windows(tumult.read(str(mapped)), LINES_SED)
            m['snaps'][snap] = dict(fe=fe, nlte=found[0].name, lte=lte[0].name, spectra=sp)
        res[model] = m
        print(model, sorted(m['spectra']), {s: (v['nlte'], v['lte']) for s, v in m['snaps'].items()})
    with open(out, 'wb') as f:
        pickle.dump(dict(lines=LINES_SED, models=res), f)


def smooth(x, y, xo):
    """y(x) convolved to R_INST and interpolated onto xo."""
    lx = np.log(x)
    step = np.median(np.diff(lx))
    g = np.arange(lx[0], lx[-1], step)
    yg = np.interp(g, lx, y)
    sig = 1 / R_INST / 2.3548 / step
    k = np.exp(-0.5 * (np.arange(-int(5 * sig) - 1, int(5 * sig) + 2) / sig) ** 2)
    yc = np.convolve(np.pad(yg, len(k) // 2, mode='edge'), k / k.sum(), mode='valid')
    return np.interp(np.log(xo), g, yc)


def read_obs(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2]


def ew_pair(xo, yo, eo, xm, ym, lam0, hw):
    """(EW_obs, EW_model) in mA over the core after the linear continuum fit."""
    sel = np.abs(xo - lam0) < WIN - 0.2
    xo, yo, eo = xo[sel], yo[sel], eo[sel]
    if len(xo) < 20 or not np.all(np.isfinite(eo)):
        return np.nan, np.nan
    mo = smooth(xm, ym, xo)
    core = np.abs(xo - lam0) <= max(hw, HW_MIN)
    w = (~core) & (eo > 0)
    if w.sum() < 10:
        return np.nan, np.nan
    A = np.vstack([mo[w], mo[w] * (xo[w] - lam0)]).T / eo[w, None]
    c, *_ = np.linalg.lstsq(A, yo[w] / eo[w], rcond=None)
    yc = yo / (c[0] + c[1] * (xo - lam0))
    return (float(np.trapezoid(1 - yc[core], xo[core])) * 1e3,
            float(np.trapezoid(1 - mo[core], xo[core])) * 1e3)


def model_spectra(m):
    """kind -> {line: (x, y)}, 3D kinds averaged over the snapshots on the first one's grid;
    3d_nlte Fe lines only from snapshots with Fe departures."""
    out = dict(m['spectra'])
    snaps = list(m['snaps'].values())

    def mean(kind, have):
        res = {}
        for lam0, (x, _) in have[0]['spectra'][kind].items():
            ys = [np.interp(x, *s['spectra'][kind][lam0]) for s in have]
            res[lam0] = (x, np.mean(ys, axis=0))
        return res

    for kind in ('3d_nlte', '3d_lte', 'map_nlte'):
        have = [s for s in snaps if kind in s['spectra']]
        if have:
            out[kind] = mean(kind, have)
    with_fe = [s for s in snaps if s['fe']]
    if with_fe:
        out['3d_nlte_fe'] = mean('3d_nlte', with_fe)
    return out


def compare(win_path, obs_dir, out, tol=(150, 0.3, 0.3)):
    win_path, obs_dir = Path(win_path), Path(obs_dir)
    W = pickle.load(open(win_path, 'rb'))
    spec = pd.read_csv(win_path.parent / 'labels.csv', dtype={'specid': str, 'set': str})
    lab = pd.read_csv(win_path.parent / 's4_labels.csv', dtype={'set': str})
    lab = lab[lab.teff.notna()]
    spec = spec[spec.snr >= 30].sort_values('snr', ascending=False).drop_duplicates('OBJECT')
    rows = []
    for model, m in W['models'].items():
        near = lab[(abs(lab.teff - m['teff']) <= tol[0]) & (abs(lab.logg - m['logg']) <= tol[1])
                   & (abs(lab.mh - m['feh']) <= tol[2])]
        stars = spec.merge(near[['OBJECT', 'teff', 'logg', 'mh']], on='OBJECT',
                           suffixes=('_crude', ''))
        ms = model_spectra(m)
        kinds = [k for k in KINDS if k in ms]
        print(f'{model}: {len(stars)} stars, spectra {kinds}', flush=True)
        for _, st in stars.iterrows():
            files = list((obs_dir / st['set'] / st['stack']).glob(f'*_{st.specid}.txt'))
            if not files:
                continue
            xo, yo, eo = read_obs(files[0])
            for el, lam0, hw in W['lines']:
                row = dict(model=model, OBJECT=st.OBJECT, set=st['set'], snr=st.snr,
                           teff=st.teff, logg=st.logg, mh=st.mh, element=el, line=lam0,
                           slope=m['slope'][lam0])
                for k in kinds:
                    src = k
                    if el == 'FeII' and k.endswith('nlte') and k.replace('_nlte', '_lte') in ms:
                        src = k.replace('_nlte', '_lte')
                    elif el == 'Fe' and k == '3d_nlte':
                        if '3d_nlte_fe' not in ms:
                            continue              # no snapshot with 3D Fe departures
                        src = '3d_nlte_fe'
                    eo_, em = ew_pair(xo, yo, eo, *ms[src][lam0], lam0, hw)
                    row[f'ew_obs_{k}'], row[f'ew_{k}'] = eo_, em
                    ok = eo_ > 0 and em > 0 and row['slope'] > 0
                    row[f'dA_{k}'] = (np.log10(eo_) - np.log10(em)) / row['slope'] if ok else np.nan
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False, float_format='%.5g')
    summary(df)


def summary(df):
    """Median dA over the stars per line, then per model and spectrum: mean Mg I, mean
    and spread of Fe I over its lines, [Mg/Fe] = <Mg I> - <Fe I>."""
    cols = [c for c in df.columns if c.startswith('dA_')]
    med = df.groupby(['model', 'element', 'line'])[cols].median()
    n = df.groupby('model').OBJECT.nunique()
    for model, g in med.groupby(level=0):
        print(f'\n{model} ({n[model]} stars): median dA per line')
        print(g.droplevel(0).round(3).to_string())
        mg, fe = g.xs('Mg', level=1), g.xs('Fe', level=1)
        print('  Mg I mean ', mg.mean().round(3).to_dict())
        print('  Fe I mean ', fe.mean().round(3).to_dict())
        print('  Fe I sd   ', fe.std().round(3).to_dict())
        print('  [Mg/Fe]   ', (mg.mean() - fe.mean()).round(3).to_dict())


if __name__ == '__main__':
    if sys.argv[1:2] == ['extract'] and len(sys.argv) == 3:
        extract(sys.argv[2])
    elif sys.argv[1:2] == ['compare'] and len(sys.argv) in (5, 8):
        tol = tuple(float(v) for v in sys.argv[5:8]) or (150, 0.3, 0.3)
        compare(sys.argv[2], sys.argv[3], sys.argv[4], tol)
    else:
        sys.exit(__doc__)
