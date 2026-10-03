"""The pilot's stage-2 SEDs per Mg line: 3D NLTE on the snapshot's own departures
against 3D on the mapped <3D> departures, 3D LTE, and 1D NLTE/LTE on the <3D> average
and the MARCS twin, as abundance-like differences from 3D NLTE.

    python pilot_compare.py OUT.csv

Runs (runs/p1sed/<model>/): <snap>_3dmg (true 3D departures), <snap>_mapmg (mapped),
marcs and av3d (1D, Mg and Fe departures).  The 3D runs carry Mg departures only, so
only Mg lines are compared.  EWs and the curve-of-growth slope are those of
delta_maps.py; the slope comes from the batch-1a run nearest to the model.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import tumult
from delta_maps import RUNS, LINES_SED, ew, slopes, nearest_b1, corr
from make_nml import PILOT

MG = [l for l in LINES_SED if l[0] == 'Mg']


def ews(folder):
    r = tumult.read(str(folder))
    return {lam0: (ew(r, lam0, hw), ew(r, lam0, hw, LTE=True)) for _, lam0, hw in MG}


def main(out):
    names = sorted(p.name for p in (RUNS / 'b1/hrs').iterdir())
    rows = []
    for folder, marcs, av, feh, alpha in PILOT:
        model = Path(folder).name.replace('_new', '')
        top = RUNS / 'p1sed' / model
        teff = 5777 if model.startswith('t5777') else int(model[1:3]) * 100
        logg = int(model.split('g')[1][:2]) / 10
        s = slopes(f'b1/hrs/{nearest_b1(teff, logg, feh, names)}', MG)
        one = {k: ews(top / k) for k in ('marcs', 'av3d') if (top / k / 'save').is_dir()}
        for true in sorted(top.glob('s*_3dmg')):
            snap = true.name.split('_')[0]
            mapped = top / f'{snap}_mapmg'
            if not (true / 'save').is_dir() or not (mapped / 'save').is_dir():
                continue
            t, m = ews(true), ews(mapped)
            for lam0 in t:
                n3, l3 = t[lam0]
                row = dict(model=model, snap=snap, teff=teff, logg=logg, feh=feh, line=lam0,
                           slope=s[lam0], ew_3d_nlte=n3, ew_3d_lte=l3, ew_map_nlte=m[lam0][0],
                           ew_map_lte=m[lam0][1])
                # dA(x) = A(x) - A(3D NLTE) for the same observed EW: the bias of x
                row['dA_map'] = corr(n3, m[lam0][0], s[lam0])
                row['dA_3d_lte'] = corr(n3, l3, s[lam0])
                for k, (n1, l1) in ((k, v[lam0]) for k, v in one.items()):
                    row[f'ew_{k}_nlte'], row[f'ew_{k}_lte'] = n1, l1
                    row[f'dA_{k}_nlte'] = corr(n3, n1, s[lam0])
                    row[f'dA_{k}_lte'] = corr(n3, l1, s[lam0])
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False, float_format='%.5g')
    cols = [c for c in df.columns if c.startswith('dA_')]
    print(df[['model', 'snap', 'line'] + cols].round(3).to_string(index=False))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
