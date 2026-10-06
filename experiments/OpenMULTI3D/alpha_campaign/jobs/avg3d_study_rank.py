"""Rank the <3D> models of the averaging study by Teff_1D - Teff_3D.

    python avg3d_study_rank.py TEFF_1D TEFF_3D

Both arguments are outputs of avg3d_study_teff.py.  Teff_3D of a model is the Teff of
the snapshots' mean flux, (<Teff^4>)^(1/4); the table lists each 1D run's difference
per model and the rms over the models, best first.
"""
import sys
from collections import defaultdict

import numpy as np


def load(path):
    out = defaultdict(dict)
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        m, r, t = line.split()
        out[m][r] = float(t)
    return out


one, three = load(sys.argv[1]), load(sys.argv[2])
models = ['t5777g44m00', 't60g40m20', 't45g25m00', 't45g25m20']
t3d = {}
print('3D:')
for m in models:
    ts = np.array(list(three[m].values()))
    # Teff of the mean flux over the snapshots
    t3d[m] = (np.mean(ts**4))**0.25
    print(f'  {m} ' + ' '.join(f'{k}={v:.1f}' for k, v in three[m].items()) +
          f'  mean-flux {t3d[m]:.1f}  spread {ts.std(ddof=1):.1f}')
names = sorted({r for m in models for r in one[m]})
rows = []
for r in names:
    d = [one[m].get(r, np.nan) - t3d[m] for m in models]
    rows.append((np.sqrt(np.nanmean(np.square(d))), r, d))
rows.sort()
print(f'{"run":24s} ' + ' '.join(f'{m:>12s}' for m in models) + '   rms')
for rms, r, d in rows:
    print(f'{r:24s} ' + ' '.join(f'{x:12.1f}' for x in d) + f' {rms:6.1f}')
