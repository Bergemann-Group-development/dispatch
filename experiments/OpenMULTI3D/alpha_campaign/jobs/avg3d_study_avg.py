"""<3D> averages of the averaging study at full resolution (tumult.avg3d).

    python avg3d_study_avg.py SNAPDIR OUT MODEL FEH ALPHA SNAP [SNAP ...]

SNAPDIR holds the model's Stagger-grid-2 files (atm3d.* and mesh.*); SNAP are snapshot
numbers.  Writes OUT/<MODEL>/<surface>_<recipe>_<depth>_n<k>.atmos for every surface,
recipe and depth scale (tau only on tau500 surfaces), with k the number of pooled
snapshots, plus the opacity table (table.npz) and a summary of the tau ranges and of
the share of columns extended above their top (summary.txt).
"""
import sys
from pathlib import Path

import numpy as np
from tumult import avg3d
from tumult.hop import read_atmosphere


def main(snapdir, out, model, feh, alpha, snaps):
    files = [next(Path(snapdir).glob(f'atm3d.*_{int(s):05d}_*')) for s in snaps]
    atm = [read_atmosphere(str(f), 'stagger2') for f in files]
    out = Path(out) / model
    out.mkdir(parents=True, exist_ok=True)
    table = avg3d.OpacityTable.for_snapshots(atm, feh=feh, alpha_enhancement=alpha)
    table.save(out / 'table.npz')
    cols = avg3d.Columns.concatenate([avg3d.columns(a, table) for a in atm])
    n = len(snaps)
    lines = [f'# {model} [Fe/H]={feh} alpha={alpha} snapshots {" ".join(snaps)}',
             '# file  x_top x_bottom  extended_at_top  <ltau500>_top <ltau500>_bottom']
    for surface in avg3d.SURFACES:
        for recipe in avg3d.RECIPES:
            av = avg3d.surface_average(cols, surface, recipe)
            for depth in (('tau', 'z') if surface == 'tau500' else ('z',)):
                name = f'{surface}_{recipe}_{depth}_n{n}'
                avg3d.write_text(out / f'{name}.atmos', av, feh, depth=depth,
                                 name=f'{model}_{name}')
                lines.append(f'{name} {av["x"][0]:.3f} {av["x"][-1]:.3f} '
                             f'{av["extrapolated"][0]:.3f} {av["ltau500"][0]:.3f} '
                             f'{av["ltau500"][-1]:.3f}')
    with open(out / 'summary.txt', 'a') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], a[2], float(a[3]), float(a[4]), a[5:])
