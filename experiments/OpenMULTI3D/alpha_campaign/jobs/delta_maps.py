"""Equivalent-width corrections of batch 2: Delta_NLTE over the 1D MARCS grid (A) and
Delta_<3D> between each Stagger <3D> average and its MARCS twin (B).

    python delta_maps.py A OUT.csv      # one row per MARCS model, element and line
    python delta_maps.py B OUT.csv      # one row per <3D> average, element and line
    python delta_maps.py A2 OUT.csv     # the same from the stage-2 line-list SEDs
    python delta_maps.py B2 OUT.csv
    python delta_maps.py V OUT.csv      # Delta_<3D> in LTE, <3D> own vmic and twin vmic
    python delta_maps.py X OUT.csv      # Delta_3D in LTE: 3D snapshots against <3D>, MARCS
    python delta_maps.py N OUT.csv      # own <3D> averages (b4) against the Stagger files (b2)
    python delta_maps.py W OUT.csv      # Delta_3D in LTE against 1D on a vmic grid (b5)

Every EW is the integral of 1 - F/F_cont over the line's window around its centre, on
the R=2e5 spectra.  The windows include blends, which are the same in both spectra of
a pair; strong lines get a wider window that still cuts their wings (a core EW).

A and B read the departure runs (b2/n1d, b2/av3d, b2/twin): the NLTE and LTE spectra
of an element come from its NLTE run and carry the model atom's lines (LINES, model
atom wavelengths).  fe607a has no Fe II bound-bound transitions and these runs exclude
VALD Fe, so they have no Fe II lines.  A2 and B2 read the stage-2 SEDs (b2/sed, from
make_nml.py batch2_sed): the full line list on the Mg and Fe departures, NLTE and LTE
from one run, with Fe II (LINES_SED, VALD wavelengths).

A difference in EW becomes an abundance-like correction through the 1D LTE curve of
growth of the batch-1a run of the MARCS model (or the nearest batch-1a model): its
[X/Fe] = +-0.3 sets give s = dlog(EW)/dA, and

    dA = -(log EW_a - log EW_b) / s

is the abundance change that would give spectrum b the EW of spectrum a: the
correction to add to an abundance from b (1D LTE, or MARCS) to get a (NLTE, or <3D>).
The b1 sets use the VALD line, not the model atom's, so s is approximate where their
gf values differ.

V and X read plain LTE runs with the full line list (LINES_SED): V the <3D> average with
its own vmic (b2/av3d/<av>/lte) and with its MARCS twin's constant vmic
(b3/av3d_vmic/<av>/lte), each against the twin (dA_av = A(<3D>) - A(MARCS)); X every
finished 3D LTE snapshot (x3d/, and the pilot's p1/<model>/s*_lte) against those three
(dA_3d_av, dA_3d_avv, dA_3d_tw = A(3D) - A(1D)).

N compares the own <3D> averages with density (make_nml.py avg3d, runs b4/) with the
Stagger MULTI-format averages of batch 2: in LTE (lte runs) and, where b4 has them, in
the stage-2 SEDs (NLTE and LTE).  dA_new_old = A(own) - A(Stagger file) at the same EW;
dA_3d_* = A(own <3D>) - A(MARCS twin), the new Delta_<3D>.  It also records the run's
log tau_500 against the file's at three temperatures (the phase-3 tau-scale check).

W reads every finished 3D LTE snapshot and its 1D LTE partners on the constant-vmic grid
of make_nml.py vmic_grid (b5/vmic/<model>/{marcs,av3d}_v<vmic>): dA_3d_<kind>_v<vmic> =
A(3D) - A(1D at that vmic) per line.  Per snapshot and 1D kind it also fits the vmic at
which the Fe I lines give one abundance (the minimum over vmic, interpolated in steps of
0.05 km/s, of the spread of dA over the Fe I lines, which removes the trend with line
strength) and records Delta_3D there (columns vfit_<kind>, dA_3d_<kind>_vfit).
"""
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import tumult

CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')
RUNS = CAMPAIGN / 'runs'
# (element, centre [A, air, model atom], half window [A])
LINES = [('Mg', 5172.681, 1.5), ('Mg', 5528.416, 0.8), ('Mg', 5711.063, 0.4),
         ('Fe', 5232.939, 1.0), ('Fe', 5269.537, 1.0), ('Fe', 5250.209, 0.3),
         ('Fe', 5383.369, 0.5), ('Fe', 6136.615, 0.5), ('Fe', 6252.555, 0.5),
         ('Fe', 6430.845, 0.4), ('Fe', 6593.869, 0.4)]
LINES_SED = [('Mg', 5172.684, 1.5), ('Mg', 5528.405, 0.8), ('Mg', 5711.088, 0.4),
             ('Fe', 5232.940, 1.0), ('Fe', 5269.537, 1.0), ('Fe', 5250.209, 0.3),
             ('Fe', 5383.369, 0.5), ('Fe', 6136.615, 0.5), ('Fe', 6252.555, 0.5),
             ('Fe', 6430.846, 0.4), ('Fe', 6593.870, 0.4),
             ('FeII', 5234.625, 0.5), ('FeII', 6247.557, 0.4), ('FeII', 6432.680, 0.4),
             ('FeII', 6456.383, 0.4), ('FeII', 6516.080, 0.4)]
NAME = re.compile(r'[ps](?P<teff>\d+)_g(?P<logg>[+-]\d\.\d)_m[\d.]+_t\d\d_z(?P<feh>[+-]\d\.\d\d)$')


def ew(run, lam0, hw, **kw):
    s = run.crop(xmin=lam0 - hw, xmax=lam0 + hw, **kw)
    x, y = np.asarray(s.xdata.value), np.asarray(s.ydata.value)
    return float(np.trapezoid(1 - y, x)) * 1e3            # mA


def set_index(run, element, da):
    for i, s in enumerate(run.sets):
        if s['offsets'] == {element: da}:
            return i
    raise KeyError(f'no {element} {da:+} set')


def slopes(b1_run, lines=LINES):
    """dlog(EW)/dA of every line from the +-0.3 sets of a batch-1a run (Fe sets for
    Fe II)."""
    r = tumult.read(str(RUNS / b1_run))
    out = {}
    for el, lam0, hw in lines:
        el = el[:2]
        hi = ew(r, lam0, hw, index=set_index(r, el, 0.3))
        lo = ew(r, lam0, hw, index=set_index(r, el, -0.3))
        out[lam0] = np.log10(hi / lo) / 0.6 if hi > 0 and lo > 0 else np.nan
    return out


def nlte_lte(folder, lines=LINES):
    """{line: (EW NLTE, EW LTE)}: of the run's element for a departure run (folder
    ending in mg or fe), of every line for a stage-2 SED."""
    r = tumult.read(str(RUNS / folder))
    sel = lines
    if lines is LINES:
        el = 'Mg' if folder.endswith('mg') else 'Fe'
        sel = [l for l in lines if l[0] == el]
    for e, lam0, hw in sel:
        y = np.asarray(r.crop(xmin=lam0 - hw, xmax=lam0 + hw).ydata.value)
        PEAK[lam0] = max(PEAK.get(lam0, 0.0), float(y.max()))
    return {lam0: (e, ew(r, lam0, hw), ew(r, lam0, hw, LTE=True)) for e, lam0, hw in sel}


# Highest normalised NLTE flux in each line window over the runs a row has read
# (reset per row by the one_* functions); rows above 1 + EMISSION_TOL are flagged.
PEAK = {}
EMISSION_TOL = 1e-3


def flag_rows(df):
    """Fe II as LTE, emission flags (phase-4 decisions 4 and 8).

    fe607a has no Fe II bound-bound transitions, and its excited Fe II super-levels run
    40-60 % ahead of the ground in metal-poor giants (phase 3), so its Fe II corrections
    are not used: in Fe II rows every NLTE EW and correction takes its LTE counterpart
    (Delta_NLTE = 0) and feii_lte is set.  emission marks lines whose NLTE spectrum
    rises above the continuum in the window (nlte_peak), the high-excitation Fe I
    lines of the [Fe/H] = -4 giants; they are kept as a model prediction."""
    if 'element' not in df:
        return df
    fe2 = df['element'] == 'FeII'
    df['feii_lte'] = fe2
    for c in [c for c in df.columns if 'nlte' in c and c != 'nlte_peak']:
        if c.startswith('dA_nlte'):
            df.loc[fe2, c] = 0.0
        elif c.replace('nlte', 'lte') in df:
            df.loc[fe2, c] = df.loc[fe2, c.replace('nlte', 'lte')]
        else:
            df.loc[fe2, c] = np.nan
    if 'nlte_peak' in df:
        df['emission'] = (df['nlte_peak'] > 1 + EMISSION_TOL) & ~fe2
    return df


def corr(ew_a, ew_b, s):
    ok = ew_a > 0 and ew_b > 0 and np.isfinite(s) and s > 0
    return -(np.log10(ew_a) - np.log10(ew_b)) / s if ok else np.nan


def labels(name):
    m = NAME.match(name)
    return float(m['teff']), float(m['logg']), float(m['feh'])


def one_a(args):
    PEAK.clear()
    name, sed = args
    lines = LINES_SED if sed else LINES
    try:
        s = slopes(f'b1/hrs/{name}', lines)
        res = (nlte_lte(f'b2/sed/{name}', lines) if sed else
               {**nlte_lte(f'b2/n1d/{name}_mg'), **nlte_lte(f'b2/n1d/{name}_fe')})
        return [dict(model=name, element=e, line=lam0, ew_lte=el_, ew_nlte=en, slope=s[lam0],
                     dA_nlte=corr(en, el_, s[lam0]), nlte_peak=PEAK.get(lam0, np.nan))
                for lam0, (e, en, el_) in res.items()]
    except Exception as exc:
        print(f'skip {name}: {type(exc).__name__}: {exc}', flush=True)
        return []


def nearest_b1(teff, logg, feh, names):
    lab = np.array([labels(n) for n in names])
    d = ((lab[:, 0] - teff) / 250) ** 2 + ((lab[:, 1] - logg) / 0.5) ** 2 + ((lab[:, 2] - feh) / 0.5) ** 2
    return names[int(np.argmin(d))]


def one_b(args):
    PEAK.clear()
    row, b1_names, sed = args
    av, teff, logg, feh, alpha, av_run, tw_lte, tw_mg, tw_fe, tw_model = row
    teff, logg, feh = float(teff), float(logg), float(feh)
    lines = LINES_SED if sed else LINES
    try:
        b1 = nearest_b1(teff, logg, feh, b1_names)
        s = slopes(f'b1/hrs/{b1}', lines)
        if sed:
            a = nlte_lte(f'b2/sed/av3d_{av}', lines)
            tw = (f'b2/sed/twin_{av}' if tw_mg.startswith('b2/twin/') else
                  f'b2/sed/{Path(tw_mg).name[:-3]}' if tw_mg.startswith('b2/n1d/') else '-')
            b = nlte_lte(tw, lines) if tw != '-' else {k: (v[0], np.nan, np.nan) for k, v in a.items()}
        else:
            a = {**nlte_lte(f'{av_run}/mg'), **nlte_lte(f'{av_run}/fe')}
            b = {}
            for t in (tw_mg, tw_fe):
                b.update(nlte_lte(t) if t != '-' else {})
            b = {k: b.get(k, (v[0], np.nan, np.nan)) for k, v in a.items()}
        rows = []
        for lam0, (e, an, al) in a.items():
            _, bn, bl = b[lam0]
            rows.append(dict(av3d=av, teff=teff, logg=logg, feh=feh, alpha=float(alpha),
                             twin=tw_model, slope_from=b1, element=e, line=lam0,
                             ew_av_lte=al, ew_av_nlte=an, ew_tw_lte=bl, ew_tw_nlte=bn,
                             slope=s[lam0], dA_3d_lte=corr(al, bl, s[lam0]),
                             dA_3d_nlte=corr(an, bn, s[lam0]),
                             dA_nlte_av=corr(an, al, s[lam0]), dA_nlte_tw=corr(bn, bl, s[lam0]),
                             nlte_peak=PEAK.get(lam0, np.nan)))
        return rows
    except Exception as exc:
        print(f'skip {av}: {type(exc).__name__}: {exc}', flush=True)
        return []


def lte_ews(folder, lines=LINES_SED):
    r = tumult.read(str(RUNS / folder))
    return {lam0: ew(r, lam0, hw, LTE=True) for _, lam0, hw in lines}


def pair_rows():
    return [l.split() for l in (CAMPAIGN / 'work/b2/pairs.txt').read_text().splitlines()
            if not l.startswith('#')]


def one_v(args):
    row, b1_names = args
    av, teff, logg, feh, alpha, av_run, tw_lte = row[:7]
    if tw_lte == '-':
        return []
    try:
        s = slopes(f'b1/hrs/{nearest_b1(float(teff), float(logg), float(feh), b1_names)}', LINES_SED)
        a, v, b = lte_ews(f'{av_run}/lte'), lte_ews(f'b3/av3d_vmic/{av}/lte'), lte_ews(tw_lte)
        return [dict(av3d=av, teff=float(teff), logg=float(logg), feh=float(feh), element=e,
                     line=lam0, ew_av=a[lam0], ew_avv=v[lam0], ew_tw=b[lam0], slope=s[lam0],
                     dA_av=corr(a[lam0], b[lam0], s[lam0]), dA_avv=corr(v[lam0], b[lam0], s[lam0]))
                for e, lam0, _ in LINES_SED]
    except Exception as exc:
        print(f'skip {av}: {type(exc).__name__}: {exc}', flush=True)
        return []


def one_x(args):
    run, row, b1_names = args
    av, teff, logg, feh, alpha, av_run, tw_lte = row[:7]
    try:
        s = slopes(f'b1/hrs/{nearest_b1(float(teff), float(logg), float(feh), b1_names)}', LINES_SED)
        x = lte_ews(run)
        one = {'av': lte_ews(f'{av_run}/lte')}
        if (RUNS / f'b3/av3d_vmic/{av}/lte/save').is_dir():
            one['avv'] = lte_ews(f'b3/av3d_vmic/{av}/lte')
        if tw_lte != '-':
            one['tw'] = lte_ews(tw_lte)
        rows = []
        for e, lam0, _ in LINES_SED:
            r = dict(run=run, av3d=av, teff=float(teff), logg=float(logg), feh=float(feh),
                     element=e, line=lam0, ew_3d=x[lam0], slope=s[lam0])
            for k, v in one.items():
                r[f'ew_{k}'] = v[lam0]
                r[f'dA_3d_{k}'] = corr(x[lam0], v[lam0], s[lam0])
            rows.append(r)
        return rows
    except Exception as exc:
        print(f'skip {run}: {type(exc).__name__}: {exc}', flush=True)
        return []


def tau_check(av):
    """Run log tau_500 of the own-average LTE run at the temperatures where the file
    has log tau = -3, -1, 0, and the run's top."""
    f = np.loadtxt(CAMPAIGN / f'work/avg3d/{av}.atmos', skiprows=3)
    a = tumult.read(str(RUNS / f'b4/av3d/{av}/lte')).atmos
    lt, t = np.ravel(a.ltau), np.ravel(a.temp)
    out = {'run_top': float(lt.min()), 'file_top': float(f[:, 0].min())}
    for x in (-3, -1, 0):
        tx = np.interp(x, f[:, 0], f[:, 1])
        k = np.argmin(np.abs(t - tx) + 1e3 * (lt < x - 2))      # the branch near x
        out[f'run_ltau_at_file{x}'] = float(lt[k])
    return out


def one_n(args):
    PEAK.clear()
    row, b1_names = args
    av, teff, logg, feh, alpha, av_run, tw_lte, tw_mg = row[:8]
    if not (RUNS / f'b4/av3d/{av}/lte/save').is_dir():
        return []
    try:
        s = slopes(f'b1/hrs/{nearest_b1(float(teff), float(logg), float(feh), b1_names)}', LINES_SED)
        new, old = lte_ews(f'b4/av3d/{av}/lte'), lte_ews(f'{av_run}/lte')
        tw = lte_ews(tw_lte) if tw_lte != '-' else None
        sed = None
        if (RUNS / f'b4/sed/av3d_{av}/save').is_dir():
            sed = (nlte_lte(f'b4/sed/av3d_{av}', LINES_SED), nlte_lte(f'b2/sed/av3d_{av}', LINES_SED))
            twn = (f'b2/sed/twin_{av}' if tw_mg.startswith('b2/twin/') else
                   f'b2/sed/{Path(tw_mg).name[:-3]}' if tw_mg.startswith('b2/n1d/') else None)
            twn = nlte_lte(twn, LINES_SED) if twn else None
        tau = tau_check(av)
        rows = []
        for e, lam0, _ in LINES_SED:
            r = dict(av3d=av, teff=float(teff), logg=float(logg), feh=float(feh), element=e,
                     line=lam0, slope=s[lam0], ew_new=new[lam0], ew_old=old[lam0],
                     dA_new_old_lte=corr(new[lam0], old[lam0], s[lam0]), **tau)
            if tw:
                r['dA_3d_lte_new'] = corr(new[lam0], tw[lam0], s[lam0])
                r['dA_3d_lte_old'] = corr(old[lam0], tw[lam0], s[lam0])
            if sed:
                (_, nn, nl), (_, on, ol) = sed[0][lam0], sed[1][lam0]
                r.update(ew_new_nlte=nn, ew_old_nlte=on, ew_new_sl=nl, nlte_peak=PEAK.get(lam0, np.nan),
                         dA_new_old_nlte=corr(nn, on, s[lam0]),
                         dA_nlte_new=corr(nn, nl, s[lam0]), dA_nlte_old=corr(on, ol, s[lam0]))
                if twn:
                    r['dA_3d_nlte_new'] = corr(nn, twn[lam0][1], s[lam0])
                    r['dA_3d_nlte_old'] = corr(on, twn[lam0][1], s[lam0])
            rows.append(r)
        return rows
    except Exception as exc:
        print(f'skip {av}: {type(exc).__name__}: {exc}', flush=True)
        return []


VMIC_GRID = (0.5, 1.0, 1.5, 2.0, 2.5)


def one_w(args):
    run, row, b1_names = args
    av, teff, logg, feh = row[:4]
    try:
        s = slopes(f'b1/hrs/{nearest_b1(float(teff), float(logg), float(feh), b1_names)}', LINES_SED)
        x = lte_ews(run)
        rows = {lam0: dict(run=run, av3d=av, teff=float(teff), logg=float(logg), feh=float(feh),
                           element=e, line=lam0, ew_3d=x[lam0], slope=s[lam0])
                for e, lam0, _ in LINES_SED}
        for kind in ('marcs', 'av3d'):
            grid = [v for v in VMIC_GRID if (RUNS / f'b5/vmic/{av}/{kind}_v{v:.1f}/save').is_dir()]
            if len(grid) < 2:
                continue
            da = {}
            for v in grid:
                one = lte_ews(f'b5/vmic/{av}/{kind}_v{v:.1f}')
                for lam0 in rows:
                    da.setdefault(lam0, []).append(corr(x[lam0], one[lam0], s[lam0]))
                    rows[lam0][f'dA_3d_{kind}_v{v:.1f}'] = da[lam0][-1]
            fe1 = np.array([da[lam0] for e, lam0, _ in LINES_SED if e == 'Fe'])
            fine = np.arange(grid[0], grid[-1] + 1e-6, 0.05)
            spread = [np.nanstd([np.interp(v, grid, d) for d in fe1]) for v in fine]
            vfit = float(fine[int(np.nanargmin(spread))])
            for lam0 in rows:
                rows[lam0][f'vfit_{kind}'] = vfit
                rows[lam0][f'dA_3d_{kind}_vfit'] = float(np.interp(vfit, grid, da[lam0]))
        return list(rows.values())
    except Exception as exc:
        print(f'skip {run}: {type(exc).__name__}: {exc}', flush=True)
        return []


def lte3d_runs():
    """(run folder, <3D> name) of every finished 3D LTE run: x3d/<av>/s*_lte and the
    pilot's p1/<model>/s*_lte."""
    out = []
    for top in ('x3d', 'p1'):
        for run in sorted((RUNS / top).glob('t*/s*_lte')):
            log = CAMPAIGN / 'logs' / f'{top}_{run.parent.name}_{run.name}.log'
            if log.exists() and 'OpenMULTI3D terminated' in log.read_text(errors='replace'):
                out.append((f'{top}/{run.parent.name}/{run.name}', run.parent.name))
    return out


def main(which, out, nproc):
    names = sorted(p.name for p in (RUNS / 'b1/hrs').iterdir())
    sed = which.endswith('2')
    with Pool(nproc) as pool:
        if which == 'W':
            pairs = {p[0]: p for p in pair_rows()}
            todo = [(run, pairs[av], names) for run, av in lte3d_runs() if av in pairs]
            res = pool.map(one_w, todo, chunksize=1)
            df = pd.DataFrame([r for rows in res for r in rows])
        elif which == 'N':
            res = pool.map(one_n, [(p, names) for p in pair_rows()], chunksize=1)
            df = pd.DataFrame([r for rows in res for r in rows])
        elif which == 'V':
            res = pool.map(one_v, [(p, names) for p in pair_rows()], chunksize=2)
            df = pd.DataFrame([r for rows in res for r in rows])
        elif which == 'X':
            pairs = {p[0]: p for p in pair_rows()}
            todo = [(run, pairs[av], names) for run, av in lte3d_runs() if av in pairs]
            res = pool.map(one_x, todo, chunksize=1)
            df = pd.DataFrame([r for rows in res for r in rows])
        elif which.startswith('A'):
            res = pool.map(one_a, [(n, sed) for n in names], chunksize=4)
            df = pd.DataFrame([r for rows in res for r in rows])
            lab = np.array([labels(n) for n in df['model']])
            df.insert(1, 'teff', lab[:, 0])
            df.insert(2, 'logg', lab[:, 1])
            df.insert(3, 'feh', lab[:, 2])
        else:
            pairs = [l.split() for l in (CAMPAIGN / 'work/b2/pairs.txt').read_text().splitlines()
                     if not l.startswith('#')]
            res = pool.map(one_b, [(p, names, sed) for p in pairs], chunksize=2)
            df = pd.DataFrame([r for rows in res for r in rows])
    df = flag_rows(df)
    df.to_csv(out, index=False, float_format='%.5g')
    print(f'{len(df)} rows -> {out}')


if __name__ == '__main__':
    import os
    if sys.argv[1:2] in (['A'], ['B'], ['A2'], ['B2'], ['V'], ['X'], ['N'], ['W']):
        main(sys.argv[1], sys.argv[2], int(os.environ.get('SLURM_CPUS_PER_TASK', 1)))
    else:
        sys.exit(__doc__)
