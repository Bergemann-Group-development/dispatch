#!/usr/bin/env python3
"""Normalised spectra of every set of the atom-validation runs, in compact files.

For each run of ``runs_<batch>.csv`` this writes ``extract/<batch>/<name>.npz`` with the
air wavelengths, the set manifest (kind, element, da) and each set's spectrum divided by
the unperturbed continuum: the flux for the benchmark runs, the intensity at every mu
for the CLV runs. Run it in the run directory (a compute node, with tumult importable).

    python extract.py BATCH
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tumult.run import Run


def extract(folder, clv):
    run = Run.read(folder)
    kinds = np.array([s['kind'] for s in run.sets])
    da = np.array([next(iter(s['offsets'].values())) if s['offsets'] else 0.0 for s in run.sets])
    out = dict(lam=run.lam, kind=kinds, da=da)
    if clv:
        cnt = run.ie_cnt                                        # (nlam, nang)
        out['mu'] = np.asarray(run.vec)[:, 2]                    # direction cosine of each ray
        out['spec'] = np.stack([run._ie_sets[..., i] / cnt for i in range(len(kinds))], 0).astype(np.float32)
        out['spec'] = np.moveaxis(out['spec'], 2, 0)            # (nmu, nset, nlam)
        out['cont'] = cnt.T.astype(np.float32)
    else:
        out['spec'] = (run._flux_sets / run.flux_cnt[:, None]).T.astype(np.float32)   # (nset, nlam)
        out['cont'] = run.flux_cnt.astype(np.float32)
    return out


def main(batch):
    runs = pd.read_csv(f'runs_{batch}.csv')
    dest = Path('extract') / batch
    dest.mkdir(parents=True, exist_ok=True)
    for r in runs.itertuples():
        folder = Path('data') / Path(r.nml).with_suffix('')
        target = dest / f'{r.name}.npz'
        try:
            if target.exists() and target.stat().st_mtime > (folder / 'sets.txt').stat().st_mtime:
                continue
            np.savez_compressed(target, **extract(folder, r.kind == 'clv'))
        except Exception as err:                                # report and go on: one failed run
            print(f'{r.name}: {type(err).__name__}: {err}')


if __name__ == '__main__':
    main(sys.argv[1])
