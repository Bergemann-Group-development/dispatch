"""Diagnostics of two phase-2 anomalies in the 1D Fe NLTE runs.

    python anomalies.py feii OUT.csv          # Fe II: matched levels and their departures
    python anomalies.py emission RUN [RUN..]  # NLTE flux above the continuum
    python anomalies.py levels RUN L [L ..]   # log b of levels L (1-based) and T per depth

feii: every Fe II line of delta_maps.LINES_SED takes, in a stage-2 SED, the departure
coefficients of the fe607a levels nearest in energy (HOP find_closest; level_match.py).
For every batch-2 MARCS model this writes, per line, the matched lower and upper levels,
their energy mismatch, and log b of both levels and of the Fe II ground level at
log tau_500 = -3, -2, -1, 0 from the model's Fe departure run (b2/n1d/<model>_fe), next to
the line's Delta_NLTE from delta_A2.csv.  For a weak line with the source function
S = b_u/b_l B, the NLTE strengthening follows log(b_u/b_l) where the line forms.

emission: the wavelengths where a run's normalised NLTE flux exceeds 1, the LTE flux
there, the model-atom line nearest to each peak (the departure runs' NLTE spectra carry
the atom's Fe I lines), and log b of its levels over depth.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import tumult
from tumult.atom import AtomClass
from delta_maps import RUNS, LINES_SED

CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')
ATOM = CAMPAIGN / 'work/input/atoms/atom.fe607a'
VALD = CAMPAIGN / 'work/input/linelists/master_linelists/vald_atoms_abo.txt'
HC = 12398.419843320026        # eV A
LTAU = (-3.0, -2.0, -1.0, 0.0)


def atom_keys():
    """fe607a levels as (key, stage, label): key = 100 x stage + energy above the stage's
    ground level in eV, as HOP's find_closest compares them."""
    a = AtomClass(str(ATOM))
    ev = np.array([l.ev for l in a.level])
    ion = np.array([l.ion for l in a.level])
    keys, ground = ev.copy(), ev[0]
    for k in range(len(ev)):
        if k and ion[k] > ion[k - 1]:
            ground = ev[k]
        keys[k] = 100 * ion[k] + ev[k] - ground
    return a, keys, ion


def vald_block(species):
    """(lambda, E_low, log gf) of a VALD species block."""
    lines = VALD.read_text().splitlines()
    i = lines.index(f"'{species}'")
    n = int(lines[i + 1].split()[1])
    return np.array([[float(x) for x in l.split()[:3]] for l in lines[i + 2:i + 2 + n]])


def departures(run):
    """log tau_500 (ascending) and b (depth, level) of a 1D run's Fe departures."""
    r = tumult.read(str(run))
    pop = np.asarray(r.read_patch_save('pop_Fe', folder=Path(run) / 'save', zdim=2))
    b = pop[0, 0, :, :, 1, 0]
    ltau = np.asarray(r.atmos.ltau, float).ravel()
    order = np.argsort(ltau)
    return ltau[order], b[order]


def at(ltau, b, level):
    return [float(np.log10(np.interp(t, ltau, b[:, level]))) for t in LTAU]


def feii(out):
    _, keys, ion = atom_keys()
    ground2 = int(np.where(ion == 2)[0][0])
    vald = vald_block('Fe II')
    match = []
    for el, lam0, _ in (l for l in LINES_SED if l[0] == 'FeII'):
        lam, elo, lgf = vald[np.abs(vald[:, 0] - lam0).argmin()]
        klo, kup = 200 + elo, 200 + elo + HC / lam
        il, iu = np.abs(keys - klo).argmin(), np.abs(keys - kup).argmin()
        match.append(dict(line=round(lam0, 1), elow=elo, eup=elo + HC / lam, loggf=lgf,
                          lo=il + 1, up=iu + 1, dE_lo=keys[il] - klo, dE_up=keys[iu] - kup))
    print(pd.DataFrame(match).round(3).to_string(index=False))
    d = pd.read_csv(RUNS / 'b2/delta_A2.csv')
    d = d[d['element'] == 'FeII']
    rows = []
    for model, g in d.groupby('model'):
        try:
            ltau, b = departures(RUNS / f'b2/n1d/{model}_fe')
        except Exception as exc:
            print(f'skip {model}: {exc}', flush=True)
            continue
        bg = at(ltau, b, ground2)
        for m in match:
            sel = g[abs(g['line'] - m['line']) < 0.1]
            if not len(sel):
                continue
            r = sel.iloc[0]
            bl, bu = at(ltau, b, m['lo'] - 1), at(ltau, b, m['up'] - 1)
            row = dict(model=model, teff=r['teff'], logg=r['logg'], feh=r['feh'],
                       line=m['line'], dE_lo=m['dE_lo'], dE_up=m['dE_up'],
                       dA_nlte=r['dA_nlte'], ew_lte=r['ew_lte'], ew_nlte=r['ew_nlte'])
            for t, x, y, z in zip(LTAU, bl, bu, bg):
                row[f'lb_lo_{t:+.0f}'], row[f'lb_up_{t:+.0f}'], row[f'lb_g2_{t:+.0f}'] = x, y, z
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False, float_format='%.5g')
    df['ratio_m1'] = df['lb_up_-1'] - df['lb_lo_-1']
    df['ratio_m2'] = df['lb_up_-2'] - df['lb_lo_-2']
    print(f'{len(df)} rows -> {out}')
    print(df.groupby(['feh', 'line'])[['dA_nlte', 'lb_lo_-1', 'lb_up_-1', 'ratio_m1',
                                        'ratio_m2', 'lb_g2_-1']].median().round(3).to_string())
    ok = df['dA_nlte'].notna()
    for c in ('ratio_m1', 'ratio_m2', 'lb_lo_-1'):
        print(f'corr(dA_nlte, {c}) = {np.corrcoef(df.loc[ok, "dA_nlte"], df.loc[ok, c])[0, 1]:.2f}')


def emission(runs):
    a, keys, ion = atom_keys()
    lines = [(t.lam0, t.i, t.j) for t in a.trans[:a.nline]]
    lam_atom = np.array([l[0] for l in lines])
    for run in runs:
        r = tumult.read(str(run))
        s = r.crop(xmin=r.lam.min(), xmax=r.lam.max())
        x, f = np.asarray(s.xdata.value), np.asarray(s.ydata.value)
        fl = np.asarray(r.crop(xmin=r.lam.min(), xmax=r.lam.max(), LTE=True).ydata.value)
        ltau, b = departures(run)
        print(f'== {run}: {np.sum(f > 1)} points above the continuum, max {f.max():.4f}')
        above = np.where(f > 1.005)[0]
        peaks = [k for k in above if f[k] == f[max(k - 3, 0):k + 4].max()]
        for k in sorted(peaks, key=lambda k: -f[k])[:8]:
            j = np.abs(lam_atom - x[k]).argmin()        # both in air
            lam0, lo, up = lines[j]
            lo, up = sorted((lo, up), key=lambda i: keys[i])
            print(f'  {x[k]:10.3f} NLTE {f[k]:.4f} LTE {fl[k]:.4f}; atom line {lam0:.3f} '
                  f'levels {lo + 1} ({keys[lo] - 100:.3f} eV) - {up + 1} ({keys[up] - 100:.3f} eV)')
            for t in (-4, -3, -2, -1, 0):
                bl = np.interp(t, ltau, b[:, lo])
                bu = np.interp(t, ltau, b[:, up])
                print(f'      log tau {t:+d}: log b_lo {np.log10(bl):+.3f}  log b_up {np.log10(bu):+.3f}'
                      f'  log(b_up/b_lo) {np.log10(bu / bl):+.3f}')


def levels(run, which):
    """log tau_500, T and log b of the given levels at every depth point."""
    r = tumult.read(str(run))
    pop = np.asarray(r.read_patch_save('pop_Fe', folder=Path(run) / 'save', zdim=2))
    b = pop[0, 0, :, :, 1, 0]
    ltau = np.asarray(r.atmos.ltau, float).ravel()
    temp = np.asarray(r.atmos.temp, float).ravel()
    print(f'== {run}')
    print('  ltau      T  ' + '  '.join(f'lb{l:>4d}' for l in which))
    for k in np.argsort(ltau):
        print(f'{ltau[k]:+6.2f} {temp[k]:6.0f}  ' + '  '.join(f'{np.log10(b[k, l - 1]):+6.3f}' for l in which))


if __name__ == '__main__':
    if sys.argv[1:2] == ['feii'] and len(sys.argv) == 3:
        feii(sys.argv[2])
    elif sys.argv[1:2] == ['emission'] and len(sys.argv) > 2:
        emission(sys.argv[2:])
    elif sys.argv[1:2] == ['levels'] and len(sys.argv) > 3:
        levels(sys.argv[2], [int(x) for x in sys.argv[3:]])
    else:
        sys.exit(__doc__)
