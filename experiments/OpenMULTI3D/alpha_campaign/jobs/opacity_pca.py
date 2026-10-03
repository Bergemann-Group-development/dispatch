"""How compressible the batch-1b LTE opacity spectra are across (T, rho, [Fe/H]).

    python opacity_pca.py OUT.npz [N_TRAIN N_TEST]     # per composition, default 512 256

Reads N_TRAIN + N_TEST points per composition from the shards of opacity_points.py
(runs/opac/feh*/), each a spectrum of the total extinction chi on the 64347-point HRS
grid, and works with y = log10(chi / chi_cont): the line part of the opacity relative to
the continuum, which takes out the smooth T and rho dependence of the continuum.  A PCA
(SVD of the centred training matrix, all compositions together) gives the explained
variance against the number of components; each held-out point is projected on the first
k components and the reconstruction error |dy| reported, over all held-out points and
over those of the line-forming regime (4000 K <= T <= 7000 K, 1e-9 <= rho <= 1e-6
g/cm3), as the median, 99th percentile and maximum over wavelengths and points.  A
log10 error of 0.01 is a 2.3 % error in the line opacity.
"""
import sys
from pathlib import Path

import numpy as np

OPAC = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/runs/opac')
KS = (1, 2, 4, 8, 16, 32, 64, 128, 256)


def load(n_train, n_test, seed=1):
    rng = np.random.default_rng(seed)
    train, test, meta = [], [], []
    for folder in sorted(OPAC.glob('feh*')):
        shards = sorted(folder.glob('shard_*.npz'))
        per = int(np.ceil((n_train + n_test) / len(shards)))
        y, tr = [], []
        for s in shards:
            d = np.load(s)
            pick = rng.choice(len(d['temp']), min(per, len(d['temp'])), replace=False)
            chi, cont = d['chi'][pick].astype(np.float64), d['chi_cont'][pick].astype(np.float64)
            y.append(np.log10(np.maximum(chi, cont) / cont).astype(np.float32))
            tr.append(np.column_stack([d['temp'][pick], d['rho'][pick],
                                       np.full(len(pick), float(d['feh']))]))
        y, tr = np.concatenate(y), np.concatenate(tr)
        order = rng.permutation(len(y))
        train.append(y[order[:n_train]])
        test.append(y[order[n_train:n_train + n_test]])
        meta.append(tr[order[n_train:n_train + n_test]])
        print(f'{folder.name}: {len(y)} points read', flush=True)
    return np.concatenate(train), np.concatenate(test), np.concatenate(meta)


def main(out, n_train=512, n_test=256):
    train, test, meta = load(n_train, n_test)
    mean = train.mean(axis=0)
    a = train - mean
    gram = (a @ a.T).astype(np.float64)             # points x points
    w, v = np.linalg.eigh(gram)
    w, v = w[::-1], v[:, ::-1]
    keep = w > w[0] * 1e-12
    basis = (v[:, keep].T.astype(np.float32) @ a) / np.sqrt(w[keep]).astype(np.float32)[:, None]
    explained = np.cumsum(w[keep]) / w[keep].sum()
    photo = ((meta[:, 0] >= 4000) & (meta[:, 0] <= 7000) &
             (meta[:, 1] >= 1e-9) & (meta[:, 1] <= 1e-6))
    t = test - mean
    rows = []
    for k in KS:
        if k > len(basis):
            break
        b = basis[:k]
        err = np.abs(t - (t @ b.T) @ b)
        for name, sel in (('all', np.ones(len(t), bool)), ('photosphere', photo)):
            e = err[sel]
            rows.append((k, name, explained[k - 1], np.median(e), np.percentile(e, 99), e.max()))
    print(f'{len(train)} training and {len(test)} held-out points ({photo.sum()} photospheric), '
          f'{train.shape[1]} wavelengths; y = log10(chi/chi_cont), range {train.min():.2f}..{train.max():.2f}')
    print(f'{"k":>4} {"set":>12} {"expl.var":>9} {"med|dy|":>9} {"p99|dy|":>9} {"max|dy|":>9}')
    for k, name, ex, med, p99, mx in rows:
        print(f'{k:4d} {name:>12} {ex:9.6f} {med:9.2e} {p99:9.2e} {mx:9.2e}')
    np.savez(out, explained=explained, ks=np.array([r[0] for r in rows]),
             sets=np.array([r[1] for r in rows]), med=np.array([r[3] for r in rows]),
             p99=np.array([r[4] for r in rows]), max=np.array([r[5] for r in rows]))


if __name__ == '__main__':
    if len(sys.argv) not in (2, 4):
        sys.exit(__doc__)
    main(sys.argv[1], *(int(x) for x in sys.argv[2:]))
