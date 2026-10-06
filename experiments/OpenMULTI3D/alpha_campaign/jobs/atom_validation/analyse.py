#!/usr/bin/env python3
"""Line-by-line abundances of the benchmark stars and the solar CLV from the atom-validation
runs, and the metrics that compare the original and the repaired atoms.

Equivalent widths, measured the same way in the observation and in the synthesis:

- Each observed spectrum is matched to the run's LTE synthesis at A_ref: a velocity shift
  (cross-correlation) and one Gaussian broadening (instrument, macroturbulence and
  rotation together) fitted over all Fe windows of the star.
- A line's integration window starts at lam0 +- 3 FWHM of the broadened synthesis
  (within the synthesised window) and is narrowed to the flux maxima either side of the
  line in that synthesis (``tumult.abundance.local_max_window``), so a neighbour's core
  is not counted. The same window is integrated in the observation (renormalised by the
  median ratio to the synthesis where the synthesis is continuum, within 2 %) and in
  every set of the ladder; blends in the window are therefore in both.
- The abundance is the inversion of the NLTE (and the LTE) ladder at the observed EW
  (``tumult._abundance_numerics.cog_invert``). The EW error combines the Cayrel photon
  term (SNR from the continuum residuals of the spectrum) and, for stars with several
  spectra, the spread between them.

A line counts as clean when ``tumult.abundance.screen`` passes (EW >= 5 mA, log(W/lam)
<= -4.8, blend fraction <= 0.3) and the inversion is in range. Lines whose departure
levels are more than 0.1 eV off (``select_lines.py``) are kept out of the metrics and
reported separately.

    python analyse.py LINES_DIR RUN_DIR SPECTRA_DIR OUT_DIR [stars|clv|metrics ...]

RUN_DIR holds runs_<batch>.csv, reference_abundances.csv and extract/<batch>/*.npz.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tumult._abundance_numerics import cog_invert
from tumult.abundance import local_max_window, screen
from tumult.physics import air2vac, vac2air
from tumult.spectrum import Spectrum

HERE = Path(__file__).parent
C_KMS = 299792.458
FWHM_GRID = np.arange(2.0, 22.01, 0.5)          # km/s
SHIFT_GRID = np.arange(-6.0, 6.01, 0.25)        # km/s
RENORM_MAX = 0.02
CONT_LEVEL = 0.995
MISMATCH = 'mismatch_'
N_BOOT = 1000
RNG = np.random.default_rng(1)


# --- synthetic spectra ------------------------------------------------------------------

def load_synth(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def segments(lam):
    """Index ranges of the gap-free pieces of a multi-window wavelength grid."""
    step = np.diff(lam)
    cut = np.where(step > 5 * np.median(step))[0] + 1
    return np.split(np.arange(lam.size), cut)


def broaden(lam, spec, fwhm_kms):
    """Gaussian convolution in velocity of ``spec`` (..., nlam), per gap-free piece."""
    out = np.array(spec, dtype=float, copy=True)
    sigma = fwhm_kms / 2.3548 / C_KMS
    for idx in segments(lam):
        x = np.log(lam[idx])
        if idx.size < 5:
            continue
        dx = np.median(np.diff(x))
        grid = np.arange(x[0], x[-1] + dx / 2, dx)
        half = int(np.ceil(4 * sigma / dx))
        kern = np.exp(-0.5 * (np.arange(-half, half + 1) * dx / sigma) ** 2)
        kern /= kern.sum()
        flat = out[..., idx].reshape(-1, idx.size)
        for row in range(flat.shape[0]):
            y = np.interp(grid, x, flat[row])
            pad = np.pad(y, half, mode='edge')
            flat[row] = np.interp(x, grid, np.convolve(pad, kern, mode='valid'))
        out[..., idx] = flat.reshape(out[..., idx].shape)
    return out


def set_index(synth, kind, da):
    hit = np.where((synth['kind'] == kind) & np.isclose(synth['da'], da))[0]
    return int(hit[0])


def ladder(synth, kind):
    """(da values, set indices) of one kind's curve of growth, the unperturbed set included."""
    sel = [(float(d), i) for i, (k, d) in enumerate(zip(synth['kind'], synth['da'])) if k == kind and d > -1]
    sel.sort()
    return np.array([d for d, _ in sel]), np.array([i for _, i in sel])


# --- observed spectra -------------------------------------------------------------------

def read_obs(path):
    data = np.loadtxt(path)
    lam, flux = data[:, 0], data[:, 1]
    keep = np.isfinite(flux) & (flux > 0)
    order = np.argsort(lam[keep])
    return lam[keep][order], flux[keep][order]


def on_windows(lam_obs, lam_syn):
    """Observed pixels inside the synthesised pieces."""
    keep = np.zeros(lam_obs.size, bool)
    for idx in segments(lam_syn):
        keep |= (lam_obs >= lam_syn[idx[0]]) & (lam_obs <= lam_syn[idx[-1]])
    return keep


def match(lam_obs, flux_obs, lam_syn, ref):
    """Velocity shift and broadening of the synthesis ``ref`` that best match the
    observation in the synthesised windows, and the continuum SNR of the observation."""
    keep = on_windows(lam_obs, lam_syn)
    x, y = lam_obs[keep], flux_obs[keep]
    depth = 1 - ref
    ccf = [np.corrcoef(1 - y, np.interp(x, lam_syn * (1 + v / C_KMS), depth))[0, 1] for v in SHIFT_GRID]
    shift = SHIFT_GRID[int(np.nanargmax(ccf))]
    xs = x / (1 + shift / C_KMS)
    chi = []
    for fwhm in FWHM_GRID:
        model = np.interp(xs, lam_syn, broaden(lam_syn, ref, fwhm))
        chi.append(np.mean((y - model) ** 2))
    fwhm = FWHM_GRID[int(np.argmin(chi))]
    model = np.interp(xs, lam_syn, broaden(lam_syn, ref, fwhm))
    cont = model > CONT_LEVEL
    resid = (y / model)[cont]
    snr = 1 / np.std(resid - np.median(resid)) if cont.sum() > 20 else 100.0
    return shift, fwhm, snr


def ew(lam, flux, lo, hi):
    """EW [mA] of 1 - flux over [lo, hi] (trapezoid on the given samples)."""
    sel = (lam >= lo) & (lam <= hi)
    return 1e3 * np.trapezoid(1 - flux[sel], lam[sel]) if sel.sum() > 2 else np.nan


def line_window(lam, ref_b, lam0, fwhm_kms, seg_bounds):
    hw = max(3 * fwhm_kms / C_KMS * lam0, 0.1)
    lo, hi = max(lam0 - hw, seg_bounds[0]), min(lam0 + hw, seg_bounds[1])
    return local_max_window(Spectrum(lam, ref_b), lam0, (lo, hi))


def measure_star(star, lines, runs, extract_dir, spectra_dir):
    """Per line, spectrum and run: observed and synthetic EWs on common windows."""
    files = star.spectra.split(';')
    fe_runs = runs[(runs.element == 'Fe') & (runs.atmos == 'av3d')]
    fe_run = (fe_runs if len(fe_runs) else runs).iloc[0]          # the most lines to match on
    fe = load_synth(extract_dir / f'{fe_run["name"]}.npz')
    ref_fe = fe['spec'][set_index(fe, 'lte', 0.0)]
    obs = {}
    for f in files:
        lam_o, flux_o = read_obs(spectra_dir / f)
        shift, fwhm, snr = match(lam_o, flux_o, fe['lam'], ref_fe)
        obs[f] = (lam_o / (1 + shift / C_KMS), flux_o, shift, fwhm, snr)
    rows = []
    for r in runs.itertuples():
        syn = load_synth(extract_dir / f'{r.name}.npz')
        lam = syn['lam']
        bounds = [(lam[i[0]], lam[i[-1]]) for i in segments(lam)]
        i_ref = set_index(syn, 'lte', 0.0)
        i_free = set_index(syn, 'lte', -5.0)
        sel = lines[(lines.element == r.element) & lines.benchmark]
        for f, (lam_o, flux_o, shift, fwhm, snr) in obs.items():
            spec_b = broaden(lam, syn['spec'], fwhm)
            for ln in sel.itertuples():
                seg = next((b for b in bounds if b[0] < ln.lam < b[1]), None)
                if seg is None or not (lam_o[0] < ln.lam < lam_o[-1]):
                    continue
                lo, hi = line_window(lam, spec_b[i_ref], ln.lam, fwhm, seg)
                near = (lam_o > seg[0]) & (lam_o < seg[1])
                cont = np.interp(lam_o[near], lam, spec_b[i_ref]) > CONT_LEVEL
                ratio = flux_o[near][cont] / np.interp(lam_o[near][cont], lam, spec_b[i_ref])
                norm = np.clip(np.median(ratio), 1 - RENORM_MAX, 1 + RENORM_MAX) if cont.sum() >= 5 else 1.0
                w_obs = ew(lam_o, flux_o / norm, lo, hi)
                dlam = np.median(np.diff(lam_o[(lam_o > lo) & (lam_o < hi)])) if w_obs == w_obs else np.nan
                sig_phot = 1e3 * 1.6 * np.sqrt(fwhm / C_KMS * ln.lam * dlam) / snr
                rows.append(dict(star=star.star, spectrum=f, element=r.element, ion=ln.ion, lam=ln.lam,
                                 evi=ln.evi, loggf=ln.loggf, atom=r.atom, atmos=r.atmos, a_ref=r.a_ref,
                                 lo=lo, hi=hi, norm=norm, fwhm=fwhm, shift=shift, snr=snr,
                                 w_obs=w_obs, sig_phot=sig_phot,
                                 w_free=ew(lam, spec_b[i_free], lo, hi),
                                 **{f'w_{k}{d:+.1f}': ew(lam, spec_b[i], lo, hi)
                                    for k in ('nlte', 'lte') for d, i in zip(*ladder(syn, k))}))
    return pd.DataFrame(rows)


def invert(df):
    """Combine the spectra of each (star, line, run) and invert the NLTE and LTE ladders."""
    wcols = [c for c in df if c.startswith('w_nlte') or c.startswith('w_lte')]
    keys = ['star', 'element', 'ion', 'lam', 'atom', 'atmos']
    out = []
    for key, g in df.groupby(keys):
        w = g.w_obs.to_numpy()
        ok = np.isfinite(w)
        if not ok.any():
            continue
        n = ok.sum()
        w_mean = w[ok].mean()
        spread = w[ok].std(ddof=1) if n > 1 else 0.0
        sig = np.sqrt((g.sig_phot[ok] ** 2).mean() / n + spread ** 2 / n + (0.02 * w_mean) ** 2 * (n == 1))
        first = g.iloc[0]
        row = dict(zip(keys, key), evi=first.evi, loggf=first.loggf, a_ref=first.a_ref, n_spec=n,
                   w_obs=w_mean, sig_w=sig, w_spread=spread)
        for kind in ('nlte', 'lte'):
            cols = sorted([c for c in wcols if c.startswith(f'w_{kind}')], key=lambda c: float(c[len(kind) + 2:]))
            da = np.array([float(c[len(kind) + 2:]) for c in cols])
            wsyn = g[cols].iloc[0].to_numpy(float)
            if w_mean > 0 and np.all(wsyn > 0):
                res = cog_invert(first.a_ref + da, wsyn, w_mean, sig, first.lam)
                row[f'a_{kind}'], row[f'sig_{kind}'] = res['a'], res['sigma_a']
                row[f'oor_{kind}'] = res['out_of_range']
                row['rew'] = res['rew']
        w0 = g[f'w_lte+0.0'].iloc[0]
        row['blend'] = g.w_free.iloc[0] / w0 if w0 > 0 else np.inf
        flags, clean = screen(w_mean, row['blend'], first.lam)
        row['flags'] = ';'.join(flags)
        row['clean'] = clean and not row.get('oor_nlte', True) and np.isfinite(row.get('a_nlte', np.nan))
        out.append(row)
    return pd.DataFrame(out)


# --- solar centre-to-limb variation ---------------------------------------------------

IAG_MU = (1.0, 0.99, 0.98, 0.97, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.2)
KPNO_MU = {1.0: 'kpn_o_mu1_intensity.txt', 0.2: 'kpn_o_mu02_intensity.txt'}


def read_iag(path, lam_air):
    """IAG intensity at one mu (vacuum wavelengths, not shifted to rest) on air
    wavelengths, restricted to the span of ``lam_air``."""
    from astropy.io import fits
    data = fits.getdata(path)
    lam_vac, flux = data[0], data[1]
    lo, hi = air2vac(np.array([lam_air.min() - 2, lam_air.max() + 2]))
    sel = (lam_vac > lo) & (lam_vac < hi) & np.isfinite(flux) & (flux > 0)
    return vac2air(lam_vac[sel]), flux[sel]


def clv_observations(spectra_dir, lam_air):
    """{(source, mu): (lam_air, flux)} of the IAG and KPNO intensity atlases."""
    obs = {}
    for mu in IAG_MU:
        obs[('IAG', mu)] = read_iag(spectra_dir / 'IAG' / f'solarspectrum_mu{mu:.2f}.fits', lam_air)
    for mu, f in KPNO_MU.items():
        lam, flux = read_obs(spectra_dir / 'KPNO' / f)
        keep = (lam > lam_air.min() - 2) & (lam < lam_air.max() + 2)
        obs[('KPNO', mu)] = (lam[keep], flux[keep])
    return obs


def measure_clv(lines, runs, extract_dir, spectra_dir):
    """Observed and synthetic EWs of the CLV lines at every observed mu."""
    fe_run = runs[(runs.element == 'Fe') & (runs.atmos == 'av3d')].iloc[0]
    fe = load_synth(extract_dir / f'{fe_run["name"]}.npz')
    obs = clv_observations(spectra_dir, np.concatenate([load_synth(extract_dir / f'{n}.npz')['lam']
                                                        for n in runs.name]))
    matched = {}
    for (src, mu), (lam_o, flux_o) in obs.items():
        k = int(np.argmin(np.abs(fe['mu'] - mu)))
        shift, fwhm, snr = match(lam_o, flux_o, fe['lam'], fe['spec'][k, set_index(fe, 'lte', 0.0)])
        matched[(src, mu)] = (lam_o / (1 + shift / C_KMS), flux_o, shift, fwhm, snr, k)
        print(f'CLV {src} mu={mu}: shift {shift:+.2f} km/s, FWHM {fwhm:.1f} km/s, SNR {snr:.0f}')
    rows = []
    for r in runs.itertuples():
        syn = load_synth(extract_dir / f'{r.name}.npz')
        lam = syn['lam']
        bounds = [(lam[i[0]], lam[i[-1]]) for i in segments(lam)]
        i_ref, i_free = set_index(syn, 'lte', 0.0), set_index(syn, 'lte', -5.0)
        sel = lines[(lines.element == r.element) & lines.clv]
        for (src, mu), (lam_o, flux_o, shift, fwhm, snr, k) in matched.items():
            if abs(syn['mu'][k] - mu) > 1e-3:
                continue
            spec_b = broaden(lam, syn['spec'][k], fwhm)
            for ln in sel.itertuples():
                seg = next((b for b in bounds if b[0] < ln.lam < b[1]), None)
                if seg is None or not (lam_o[0] < ln.lam < lam_o[-1]):
                    continue
                lo, hi = line_window(lam, spec_b[i_ref], ln.lam, fwhm, seg)
                near = (lam_o > seg[0]) & (lam_o < seg[1])
                cont = np.interp(lam_o[near], lam, spec_b[i_ref]) > CONT_LEVEL
                ratio = flux_o[near][cont] / np.interp(lam_o[near][cont], lam, spec_b[i_ref])
                norm = np.clip(np.median(ratio), 1 - RENORM_MAX, 1 + RENORM_MAX) if cont.sum() >= 5 else 1.0
                w_obs = ew(lam_o, flux_o / norm, lo, hi)
                pix = np.diff(lam_o[(lam_o > lo) & (lam_o < hi)])
                dlam = np.median(pix) if pix.size else np.nan
                rows.append(dict(source=src, mu=mu, element=r.element, ion=ln.ion, lam=ln.lam, evi=ln.evi,
                                 loggf=ln.loggf, atom=r.atom, atmos=r.atmos, a_ref=r.a_ref, lo=lo, hi=hi,
                                 norm=norm, fwhm=fwhm, snr=snr, w_obs=w_obs,
                                 sig_phot=1e3 * 1.6 * np.sqrt(fwhm / C_KMS * ln.lam * dlam) / snr,
                                 w_free=ew(lam, spec_b[i_free], lo, hi),
                                 **{f'w_{kd}{d:+.1f}': ew(lam, spec_b[i], lo, hi)
                                    for kd in ('nlte', 'lte') for d, i in zip(*ladder(syn, kd))}))
    return pd.DataFrame(rows)


def read_takeda(path):
    """{line code: DataFrame(mu, w)} from Takeda & UeNo's all_clvdata (W in mA)."""
    out, code, rows = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('=========='):
            if code:
                out[code] = pd.DataFrame(rows, columns=['mu', 'w'])
            code, rows = line.split('file name:')[1].split('.dat')[0].strip(), []
        elif line[:2].isdigit():
            v = line.split()
            if float(v[4]) > 0:
                rows.append((float(v[2]), float(v[4])))
    if code:
        out[code] = pd.DataFrame(rows, columns=['mu', 'w'])
    return out


def takeda_abundances(lines, runs, extract_dir, takeda_file):
    """A(mu) of every Takeda line: their EW of the line itself against the line's own
    synthetic EW (all lines minus the element removed, over the whole synthesised piece),
    interpolated linearly in mu between the run's angles."""
    data = read_takeda(takeda_file)
    rows = []
    for r in runs.itertuples():
        syn = load_synth(extract_dir / f'{r.name}.npz')
        lam = syn['lam']
        bounds = [(lam[i[0]], lam[i[-1]]) for i in segments(lam)]
        i_free = set_index(syn, 'lte', -5.0)
        order = np.argsort(syn['mu'])
        for ln in lines[(lines.element == r.element) & lines.takeda.notna()].itertuples():
            obs = data.get(ln.takeda)
            seg = next((b for b in bounds if b[0] < ln.lam < b[1]), None)
            if obs is None or seg is None or obs.empty:
                continue
            sel = (lam >= seg[0]) & (lam <= seg[1])
            for kind in ('nlte', 'lte'):
                da, idx = ladder(syn, kind)
                w_mu = np.array([[1e3 * np.trapezoid(syn['spec'][k, i_free][sel] - syn['spec'][k, i][sel], lam[sel])
                                  for i in idx] for k in range(len(syn['mu']))])           # (nmu, nda)
                for o in obs.itertuples():
                    wsyn = np.array([np.interp(o.mu, syn['mu'][order], w_mu[order, j]) for j in range(len(da))])
                    if not np.all(wsyn > 0):
                        continue
                    res = cog_invert(r.a_ref + da, wsyn, o.w, 0.0, ln.lam)
                    rows.append(dict(source='Takeda', element=r.element, ion=ln.ion, lam=ln.lam, evi=ln.evi,
                                     atom=r.atom, atmos=r.atmos, kind=kind, mu=o.mu, w_obs=o.w, a=res['a'],
                                     rew=res['rew'], oor=res['out_of_range']))
    return pd.DataFrame(rows)


# --- metrics ----------------------------------------------------------------------------

def slope(x, y, w):
    if len(x) < 4 or np.ptp(x) == 0:
        return np.nan, np.nan
    coef, cov = np.polyfit(x, y, 1, w=np.sqrt(w), cov='unscaled')
    resid = y - np.polyval(coef, x)
    scale = np.sqrt(np.sum(w * resid ** 2) / max(len(x) - 2, 1))
    return coef[0], np.sqrt(cov[0, 0]) * scale


def metrics(g):
    """Metrics of one (star, element, atom, atmos) line set (clean lines only)."""
    out = {}
    w = 1 / np.maximum(g.sig_nlte.fillna(0.3), 0.03) ** 2
    for ion, s in g.groupby('ion'):
        tag = 'I' if ion == 1 else 'II'
        ww = w[s.index]
        out[f'n_{tag}'] = len(s)
        out[f'mean_{tag}'] = np.average(s.a_nlte, weights=ww)
        out[f'sd_{tag}'] = s.a_nlte.std(ddof=1) if len(s) > 1 else np.nan
        out[f'err_{tag}'] = out[f'sd_{tag}'] / np.sqrt(len(s)) if len(s) > 1 else s.sig_nlte.iloc[0]
        out[f'lte_{tag}'] = s.a_lte.mean()
        out[f'chi_slope_{tag}'], out[f'chi_slope_err_{tag}'] = slope(s.evi.to_numpy(), s.a_nlte.to_numpy(), ww.to_numpy())
        out[f'rew_slope_{tag}'], out[f'rew_slope_err_{tag}'] = slope(s.rew.to_numpy(), s.a_nlte.to_numpy(), ww.to_numpy())
        out[f'lam_slope_{tag}'], out[f'lam_slope_err_{tag}'] = slope(s.lam.to_numpy() / 1e3, s.a_nlte.to_numpy(), ww.to_numpy())
    if 'mean_I' in out and 'mean_II' in out:
        out['ion_bal'] = out['mean_II'] - out['mean_I']
        out['ion_bal_err'] = np.hypot(out['err_I'], out['err_II'])
    return pd.Series(out)


def paired_bootstrap(a_lines, b_lines, stat):
    """Bootstrap the difference stat(b) - stat(a) over lines common to both atoms."""
    common = a_lines.index.intersection(b_lines.index)
    a, b = a_lines.loc[common], b_lines.loc[common]
    draws = []
    for _ in range(N_BOOT):
        pick = RNG.integers(0, len(common), len(common))
        draws.append(stat(b.iloc[pick]) - stat(a.iloc[pick]))
    draws = np.array(draws, float)
    return stat(b) - stat(a), np.nanstd(draws)


def main(lines_dir, run_dir, spectra_dir, out_dir, *steps):
    lines_dir, run_dir, spectra_dir, out_dir = map(Path, (lines_dir, run_dir, spectra_dir, out_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = pd.read_csv(lines_dir / 'lines.csv')
    stars = pd.read_csv(HERE / 'stars.csv', comment='#')
    runs = pd.read_csv(run_dir / 'runs_b.csv')
    steps = steps or ('stars', 'metrics')
    if 'stars' in steps:
        for s in stars.itertuples():
            target = out_dir / f'ew_{s.star}.csv'
            if target.exists():
                continue
            sel = runs[runs.star == s.star]
            have = [(run_dir / 'extract' / 'b' / f'{n}.npz').exists() for n in sel.name]
            if not all(have):
                print(f'{s.star}: {sum(have)}/{len(have)} runs extracted, skipped')
                continue
            measure_star(s, lines, sel, run_dir / 'extract' / 'b', spectra_dir).to_csv(target, index=False)
            print(f'{s.star}: measured')
    if 'metrics' in steps:
        ews = pd.concat([pd.read_csv(p) for p in sorted(out_dir.glob('ew_*.csv'))])
        res = invert(ews)
        mism = lines.set_index(['element', 'ion', 'lam'])
        flags = []
        for r in res.itertuples():
            col = MISMATCH + r.atom
            key = (r.element, r.ion, round(r.lam, 4))
            flags.append(bool(mism.loc[key, col]) if key in mism.index and col in mism else False)
        res['mismatch'] = flags
        res.to_csv(out_dir / 'line_abundances.csv', index=False)
        use = res[res.clean & ~res.mismatch]
        met = use.groupby(['star', 'element', 'atom', 'atmos']).apply(metrics, include_groups=False)
        met.to_csv(out_dir / 'metrics.csv')
        print(f'{len(res)} line results, {len(use)} clean; metrics in {out_dir / "metrics.csv"}')


if __name__ == '__main__':
    main(*sys.argv[1:])
