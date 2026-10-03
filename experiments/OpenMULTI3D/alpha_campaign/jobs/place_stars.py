"""Crude 1D LTE labels for the 4MOST HRS spectra, to place the observed stars against
the 3D Stagger grid.  These are not science labels: the grid is the batch-1a MARCS grid
(standard alpha, fixed vmic of 1 or 2 km/s, no rotation or macroturbulence), the fit
is a chi2 over the grid nodes, and only Teff, log g and [Fe/H] are free.

    python place_stars.py grid OUT.npz          # batch-1a HRS references on the observed grid
    python place_stars.py fit GRID.npz OUT.csv  # fit every observed spectrum, Gaia check

grid: the unperturbed LTE spectrum of every batch-1a HRS run (runs/b1/hrs), normalised
by its continuum, convolved with a Gaussian of FWHM lambda/R_INST in ln(lambda) per arm,
and interpolated onto the observed wavelengths (air, as in dispatch.x output).

fit: per spectrum and arm, the model times a straight line (the continuum placement of
the observed normalisation) is fitted by weighted least squares to the observed flux;
chi2 is summed over the arms.  The errors get a floor of ERR_FLOOR for model
systematics, and the Ca II H&K and H-alpha cores are masked (chromospheres, 1D LTE).
Labels are the best node and a weighted mean over all nodes with weights
exp(-(chi2 - chi2_min) / (2 s)), s = chi2_min / dof, with its spread; chi2 differences
between nodes are large, so this is close to the best node.  A free fit puts giants at
dwarf gravities (the log g - [Fe/H] - continuum degeneracy at low S/N), so a second pass
(teff_c, logg_c, feh_c) restricts log g to the photometric or asteroseismic value.

The Gaia labels: Teff from the dereddened BP-RP colour (Mucciarelli, Bellazzini & Massari
2021, the giant relation for M_G < 3, at the free-fit [Fe/H]), with E(BP-RP) = 0.432 A0;
log g from M_G, the Andrae et al. (2018) BC_G and 1 Msun, for parallax S/N > 5; and log g
from numax (scaling relation) where the 0409 photometry gives it.
"""
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')
OBS = CAMPAIGN / 'obs'
R_INST = 20000          # 4MOST HRS resolving power, all arms
ERR_FLOOR = 0.02
LOGG_WINDOW = 0.5     # second pass: grid nodes within this of the photometric log g
MASK = [(3926.0, 3980.0), (6553.0, 6573.0)]     # Ca II K and H, H-alpha core
NAME = re.compile(r'(?P<geo>[ps])(?P<teff>\d+)_g(?P<logg>[+-]\d\.\d)_m[\d.]+_t\d\d_z(?P<feh>[+-]\d\.\d\d)$')


def obs_files():
    return sorted(OBS.glob('0*/stack*/sapp_input_spectrum_*.txt'))


def read_obs(path):
    d = pd.read_csv(path, sep=r'\s+', comment='#', header=None).to_numpy()
    return d[:, 0], d[:, 1], d[:, 2]


def arms(lam):
    """Slices of the contiguous arms of a wavelength array."""
    cut = np.r_[0, np.where(np.diff(lam) > 1.0)[0] + 1, lam.size]
    return [slice(a, b) for a, b in zip(cut[:-1], cut[1:])]


def broaden_to(lam_obs, run):
    import tumult
    r = tumult.read(str(run))
    s = r.crop(xmin=float(r.lam.min()), xmax=float(r.lam.max()))
    lam, f = np.asarray(s.xdata.value), np.asarray(s.ydata.value)
    order = np.argsort(lam)
    lam, f = lam[order], f[order]
    out = np.empty(lam_obs.size, np.float32)
    for sl in arms(lam_obs):
        lo, hi = lam_obs[sl][0] - 3, lam_obs[sl][-1] + 3
        x = np.log(lam[(lam > lo) & (lam < hi)])
        step = 1 / 4e5
        u = np.arange(x[0], x[-1], step)
        fu = np.interp(u, x, f[(lam > lo) & (lam < hi)])
        fu = gaussian_filter1d(fu, 1 / R_INST / 2.3548 / step, mode='nearest')
        out[sl] = np.interp(np.log(lam_obs[sl]), u, fu)
    return out


def _grid_one(args):
    lam_obs, run = args
    try:
        return run.name, broaden_to(lam_obs, run)
    except Exception as exc:                       # a broken run is left out, and said so
        print(f'skip {run.name}: {type(exc).__name__}: {exc}', flush=True)
        return run.name, None


def grid(out, nproc):
    lam_obs = read_obs(obs_files()[0])[0]
    runs = sorted((CAMPAIGN / 'runs/b1/hrs').iterdir())
    with Pool(nproc) as pool:
        res = [r for r in pool.map(_grid_one, [(lam_obs, r) for r in runs], chunksize=4)
               if r[1] is not None]
    names = [n for n, _ in res]
    labels = np.array([[float(NAME.match(n)[k]) for k in ('teff', 'logg', 'feh')] for n in names])
    np.savez(out, lam=lam_obs, flux=np.array([f for _, f in res]), labels=labels,
             names=np.array(names))
    print(f'{len(names)} of {len(runs)} models on {lam_obs.size} pixels -> {out}')


def chi2_all(G, lam, f, e):
    """chi2 of every model row of G against one spectrum, each model times a straight
    line per arm fitted by weighted least squares; and the number of pixels used."""
    good = np.isfinite(f) & np.isfinite(e) & (e > 0) & (f > 0) & (f < 1.5)
    for lo, hi in MASK:
        good &= ~((lam > lo) & (lam < hi))
    chi2, npix = np.zeros(G.shape[0]), 0
    for sl in arms(lam):
        m = good[sl]
        if m.sum() < 100:
            continue
        g = G[:, sl][:, m].astype(np.float64)
        o = f[sl][m]
        w = 1 / (e[sl][m] ** 2 + ERR_FLOOR ** 2)
        x = lam[sl][m]
        x = (x - x.mean()) / (x.max() - x.min())
        g2 = g * g
        a11, a12, a22 = g2 @ w, g2 @ (w * x), g2 @ (w * x * x)
        b1, b2 = g @ (w * o), g @ (w * o * x)
        det = a11 * a22 - a12 * a12
        c0 = (a22 * b1 - a12 * b2) / det
        c1 = (a11 * b2 - a12 * b1) / det
        chi2 += np.sum(w * o * o) - (c0 * b1 + c1 * b2)
        npix += m.sum()
    return chi2, npix


_G = {}


def best(chi2, npix, keep):
    """Best node and weighted mean labels over the nodes in keep."""
    labels = _G['labels']
    c = np.where(keep, chi2, np.inf)
    i = int(np.argmin(c))
    s = max(c[i] / max(npix - 6, 1), 1e-12)
    w = np.exp(-(c - c[i]) / (2 * s))
    w /= w.sum()
    mean = w @ labels
    return i, mean, np.sqrt(w @ (labels - mean) ** 2)


def _fit_one(args):
    """Labels of one spectrum: free, and with log g within LOGG_WINDOW of logg_c (the
    photometric log g) when that is given."""
    path, logg_c = args
    G, lam_g, labels = _G['flux'], _G['lam'], _G['labels']
    try:
        lam, f, e = read_obs(path)
    except (ValueError, pd.errors.EmptyDataError) as exc:    # one 0408 file is empty
        print(f'skip {path.name}: {type(exc).__name__}: {exc}', flush=True)
        return None
    if lam.size != lam_g.size or not np.allclose(lam, lam_g, atol=1e-6):
        f, e = np.interp(lam_g, lam, f), np.interp(lam_g, lam, e)
    chi2, npix = chi2_all(G, lam_g, f, e)
    i, mean, sd = best(chi2, npix, np.ones(len(chi2), bool))
    ok = np.isfinite(f) & (e > 0)
    out = dict(set=path.parts[-3], stack=path.parts[-2], specid=path.stem.split('_')[-1],
               snr=float(np.nanmedian(f[ok] / e[ok])), npix=npix,
               chi2r=chi2[i] / max(npix - 6, 1),
               teff_node=labels[i, 0], logg_node=labels[i, 1], feh_node=labels[i, 2],
               teff=mean[0], logg=mean[1], feh=mean[2],
               e_teff=sd[0], e_logg=sd[1], e_feh=sd[2], model=_G['names'][i])
    if np.isfinite(logg_c):
        i, mean, _ = best(chi2, npix, np.abs(labels[:, 1] - logg_c) <= LOGG_WINDOW)
        out.update(teff_c=mean[0], logg_c=mean[1], feh_c=mean[2],
                   chi2r_c=chi2[i] / max(npix - 6, 1))
    return out


def mbm21_teff(c, feh, giant):
    """Mucciarelli, Bellazzini & Massari (2021) (BP-RP)_0 relation; NaN outside its
    colour range."""
    b = np.where(giant[:, None], [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009],
                 [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395])
    theta = b[:, 0] + b[:, 1] * c + b[:, 2] * c ** 2 + b[:, 3] * feh + b[:, 4] * feh ** 2 \
        + b[:, 5] * feh * c
    lo, hi = np.where(giant, 0.33, 0.39), np.where(giant, 1.81, 1.50)
    return np.where((c >= lo) & (c <= hi), 5040 / theta, np.nan)


def bc_g(teff):
    """Andrae et al. (2018) bolometric correction in G, 4000-8000 K."""
    dt = teff - 5772.0
    return 6.000e-02 + 6.731e-05 * dt - 6.647e-08 * dt ** 2 + 2.859e-11 * dt ** 3 - 7.197e-15 * dt ** 4


def photometry():
    """One row per spectrum id with the star's Gaia data, from every *_photometry.csv."""
    rows = []
    for csv in sorted(OBS.glob('0*_photometry.csv')):
        df = pd.read_csv(csv, dtype={c: str for c in pd.read_csv(csv, nrows=0).columns
                                     if c.startswith('PSPECUID')})
        for col in [c for c in df.columns if c.startswith('PSPECUID')]:
            part = df[df[col].notna()].copy()
            part['specid'] = part[col].str.strip()
            rows.append(part.reindex(columns=['specid', 'OBJECT', 'Plx', 'e_Plx', 'Gmag',
                                              'BPmag', 'RPmag', 'A0', 'AG', 'rgeo', 'numax']))
    # a star in several programmes keeps the row that has numax, if any
    return pd.concat(rows).sort_values('numax').drop_duplicates('specid')


def gaia_labels(res, feh):
    """Teff from the dereddened colour at [Fe/H] feh, the giant relation for M_G < 3;
    log g from M_G (parallax S/N > 5) and from numax where given."""
    a0 = res['A0'].fillna(0.0)
    c0 = res['BPmag'] - res['RPmag'] - 0.432 * a0
    dist = res['rgeo'].where(res['rgeo'] > 0, 1000 / res['Plx'])
    mg = res['Gmag'] - res['AG'].fillna(0.0) - 5 * np.log10(dist / 10)
    giant = (mg < 3.0).to_numpy()
    teff = mbm21_teff(c0.to_numpy(), np.clip(feh, -4, 0.5), giant)
    good = res['Plx'] / res['e_Plx'] > 5
    logg = np.where(good, 4.438 + 4 * np.log10(teff / 5772) + 0.4 * (mg + bc_g(teff) - 4.74), np.nan)
    # asteroseismic log g where 0409 gives numax (numax_sun = 3090 muHz)
    seis = np.log10(res['numax'] / 3090) + 4.438 + 0.5 * np.log10(teff / 5772)
    return teff, logg, seis


def fit(grid_file, out, nproc):
    """Two passes: the free fit gives [Fe/H] for the colour relation; the second pass
    restricts log g to the photometric (or asteroseismic) value."""
    g = np.load(grid_file)
    _G.update(flux=g['flux'], lam=g['lam'], labels=g['labels'], names=g['names'])
    phot = photometry().set_index('specid')
    files = obs_files()
    ids = [f.stem.split('_')[-1] for f in files]
    with Pool(nproc) as pool:
        res = pd.DataFrame([r for r in pool.map(_fit_one, [(f, np.nan) for f in files],
                                                chunksize=8) if r])
        res = res.merge(phot, left_on='specid', right_index=True, how='left')
        teff, logg, seis = gaia_labels(res, res['feh'].to_numpy())
        res['teff_phot'], res['logg_phot'], res['logg_seis'] = teff, logg, seis
        prior = dict(zip(res['specid'], np.where(np.isfinite(seis), seis, logg)))
        second = pd.DataFrame([r for r in pool.map(
            _fit_one, [(f, prior.get(i, np.nan)) for f, i in zip(files, ids)], chunksize=8) if r])
    cols = ['teff_c', 'logg_c', 'feh_c', 'chi2r_c']
    res = res.merge(second[['specid'] + [c for c in cols if c in second]], on='specid', how='left')
    res['no_A0'] = res['A0'].isna()
    res.to_csv(out, index=False, float_format='%.4g')
    print(f'{len(res)} spectra, {res["OBJECT"].nunique()} stars, '
          f'{res["teff_c"].notna().sum()} with a photometric log g -> {out}')


if __name__ == '__main__':
    nproc = int(__import__('os').environ.get('SLURM_CPUS_PER_TASK', 1))
    if sys.argv[1:2] == ['grid']:
        grid(sys.argv[2], nproc)
    elif sys.argv[1:2] == ['fit']:
        fit(sys.argv[2], sys.argv[3], nproc)
    else:
        sys.exit(__doc__)
