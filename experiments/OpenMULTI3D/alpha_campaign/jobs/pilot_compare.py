"""The pilot's stage-2 SEDs per line: 3D NLTE on the snapshot's own departures against
3D LTE, 3D on the mapped <3D> departures, and 1D NLTE/LTE on the own <3D> average and
the MARCS twin, as abundance corrections.

    python pilot_compare.py OUT.csv

3D runs (runs/p1sed/<model>/): <snap>_3dmg (true 3D Mg departures) or <snap>_3dmgfe (Mg
and Fe), with _3dmg4 / _3dmg4fe preferred (the metal-poor models' 3D Mg rerun with the
phase-4 setup), <snap>_map4mg / _map4mgfe (the own <3D> departures mapped in log tau_500).
In 1D every SED's LTE spectrum (crop(LTE=True)) is its own LTE partner.  The 3D SEDs'
LTE spectra lack the lines of the NLTE species (Mg I 5711: 2.7 against 106 mA in the Sun,
reports/2026-10-05_phase5.md, section 5), so 3D LTE comes from the plain 3D LTE run of
the snapshot (runs/p1/<model>/<snap>_lte1 with one output scattering iteration as the
SEDs, else <snap>_lte with none; VALD lines, not the model atom's).  1D: the own <3D> average's SED
(b4/sed/av3d_<model>, phase-4 departure setup) and the MARCS twin's (b6/sed/<twin> at
[Fe/H] <= -1, else p1sed/<model>/marcs).  Lines: Mg I, and Fe I where the 3D run has
Fe departures; Fe II stays LTE (phase-4 decision 4) and is left out.

Corrections, with corr(a, b) = A(a) - A(b) at the same EW (delta_maps.py):
    D_nlte_<1d>   = A(1D NLTE) - A(1D LTE)
    D_3d_lte_<1d> = A(3D LTE) - A(1D LTE)
    D_total_<1d>  = A(3D NLTE) - A(1D LTE)
    D_int_<1d>    = D_total - D_nlte - D_3d_lte  (the interaction term)
    D_nlte_3d     = A(3D NLTE) - A(3D LTE)
    D_map         = A(3D NLTE) - A(3D on mapped <3D> departures)
<1d> is av (own <3D>) or tw (MARCS twin).  The slope comes from the batch-1a run nearest
to the model.
"""
import sys
from pathlib import Path

import pandas as pd

import tumult
from delta_maps import RUNS, LINES_SED, ew, slopes, nearest_b1, corr
from make_nml import PILOT, run_name

LINES = [l for l in LINES_SED if l[0] in ('Mg', 'Fe')]


def ews(folder, lines):
    r = tumult.read(str(folder))
    return {lam0: (ew(r, lam0, hw), ew(r, lam0, hw, LTE=True)) for _, lam0, hw in lines}


def done(folder):
    return (folder / 'save').is_dir() and any((folder / 'save').glob('flux_*.bin'))


def main(out):
    names = sorted(p.name for p in (RUNS / 'b1/hrs').iterdir())
    rows = []
    for folder, marcs, av, feh, alpha in PILOT:
        model = Path(folder).name.replace('_new', '')
        top = RUNS / 'p1sed' / model
        teff = 5777 if model.startswith('t5777') else int(model[1:3]) * 100
        logg = int(model.split('g')[1][:2]) / 10
        s = slopes(f'b1/hrs/{nearest_b1(teff, logg, feh, names)}', LINES)
        tw = RUNS / f'b6/sed/{run_name(Path(marcs))}' if feh <= -1 else top / 'marcs'
        one = {k: ews(f, LINES) for k, f in (('av', RUNS / f'b4/sed/av3d_{model}'), ('tw', tw))
               if done(f)}
        for snap in sorted({p.name.split('_')[0] for p in top.glob('s*_3dmg*')}):
            # the phase-4-setup 3D Mg reruns (_3dmg4*) take precedence
            found = [f for f in (top / f'{snap}_{k}' for k in ('3dmg4fe', '3dmgfe', '3dmg4', '3dmg'))
                     if done(f)]
            if not found:
                continue
            true = found[0]
            fe = true.name.endswith('fe')
            lines = [l for l in LINES if fe or l[0] == 'Mg']
            mapped = top / f'{snap}_map4mgfe' if fe else top / f'{snap}_map4mg'
            t = ews(true, lines)
            lte3 = RUNS / f'p1/{model}/{snap}_lte1'           # one scattering iteration
            if not done(lte3):
                lte3 = RUNS / f'p1/{model}/{snap}_lte'
            r3 = tumult.read(str(lte3))
            t = {k: (v[0], ew(r3, k, hw)) for (_, k, hw), v in zip(lines, t.values())}
            m = ews(mapped, lines) if done(mapped) else {}
            for e, lam0, _ in lines:
                n3, l3 = t[lam0]
                sl = s[lam0]
                row = dict(model=model, snap=snap, teff=teff, logg=logg, feh=feh, element=e,
                           line=lam0, slope=sl, fe_3d=fe, ew_3d_nlte=n3, ew_3d_lte=l3,
                           D_nlte_3d=corr(n3, l3, sl))
                if lam0 in m:
                    row['ew_map_nlte'] = m[lam0][0]
                    row['D_map'] = corr(n3, m[lam0][0], sl)
                for k, v in one.items():
                    n1, l1 = v[lam0]
                    row[f'ew_{k}_nlte'], row[f'ew_{k}_lte'] = n1, l1
                    row[f'D_nlte_{k}'] = corr(n1, l1, sl)
                    row[f'D_3d_lte_{k}'] = corr(l3, l1, sl)
                    row[f'D_total_{k}'] = corr(n3, l1, sl)
                    row[f'D_int_{k}'] = row[f'D_total_{k}'] - row[f'D_nlte_{k}'] - row[f'D_3d_lte_{k}']
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False, float_format='%.5g')
    cols = [c for c in df.columns if c.startswith('D_')]
    print(df[['model', 'snap', 'line'] + cols].round(3).to_string(index=False))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
