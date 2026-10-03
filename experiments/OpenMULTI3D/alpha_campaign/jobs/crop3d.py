"""Cut the top of a Stagger snapshot (stagger2 format) at a horizontally averaged log tau_500.

    python crop3d.py SNAP RUN OUT [LTAU]

SNAP is an atm3d.* file (its mesh.* partner beside it), RUN a dispatch.x run of SNAP
that saved its atmosphere (an atmosphere-only make_nml.py avg3d run, or any 3D run of
the snapshot), OUT the folder for the cropped atm3d/mesh pair; LTAU defaults to -6.

The cut height is where the horizontal mean of log tau_500 over RUN's planes reaches
LTAU.  The planes above it are dropped and a plane interpolated at the cut height
(linear in height; log T, log n_e and log rho, linear velocities) becomes the new top,
so the box ends at the cut.  The bottom and the horizontal mesh stay as they are, so
the cropped and the full box share their columns and the heights below the cut.

Why: in the top planes of the 3D boxes (log tau_500 < -6, T 2000-5700 K side by side,
n_e two orders of magnitude below any 1D model) the Fe statistical equilibrium is
radiatively controlled and nearly collisionless, and the NLTE iteration there stalls
for many iterations (reports/2026-10-04_phase4.md, section 2).  The HRS lines form
below the cut (reports/2026-10-04_phase5.md, section 2 has the test).

Appends 'SNAP OUT_FILE h_cut_cm H_cm ltau_top' to OUT/crop.txt.
"""
import sys
from pathlib import Path

import numpy as np

import tumult


def read_mesh(path):
    tok = path.read_text().split()
    out, p = [], 0
    for _ in range(3):
        n = int(tok[p])
        out.append(np.array(tok[p + 1:p + 1 + n], float))
        p += 1 + n
    return out


def write_mesh(path, axes):
    with path.open('w') as f:
        for a in axes:
            f.write(f'{len(a):12d}\n')
            for i in range(0, len(a), 5):
                f.write(''.join(f'{v:16.7e}' for v in a[i:i + 5]) + '\n')


def mean_ltau_profile(run):
    """(height above the box bottom in the run's length unit, <log tau_500> per plane),
    height ascending."""
    a = tumult.read(str(run)).atmos
    zz = np.asarray(a.zz, float)
    lt = np.mean(np.asarray(a.ltau, float), axis=(0, 1))
    order = np.argsort(zz)
    return zz[order], lt[order]


def main(snap, run, out, ltau=-6.0):
    snap, out = Path(snap), Path(out)
    mesh = snap.with_name(snap.name.replace('atm3d.', 'mesh.', 1))
    x, y, z = read_mesh(mesh)
    nx, ny, nz = len(x), len(y), len(z)
    raw = np.fromfile(snap, np.float32)
    assert raw.size == 6 * nx * ny * nz, 'atm3d file holds six (nx, ny, nz) arrays'
    ne, temp, vx, vy, vz, rho = (raw[i * nx * ny * nz:(i + 1) * nx * ny * nz].reshape((nx, ny, nz), order='F')
                                 for i in range(6))
    # the top is the end with the lower density (as read_atmos_mod thick2thin decides)
    top_first = rho[0, 0, 0] < rho[0, 0, -1]
    bottom = z[-1] if top_first else z[0]
    h = np.abs(z - bottom)                          # height above the bottom plane, cm
    H = h.max()

    zr, lt = mean_ltau_profile(run)
    scale = H / zr[-1]                              # the run's top cell sits at the box top
    if not (np.all(np.diff(lt) < 0) and lt[-1] < ltau < lt[0]):
        sys.exit(f'{snap.name}: <log tau_500> not monotonic or {ltau} outside {lt[-1]:.2f}..{lt[0]:.2f}')
    h_cut = np.interp(ltau, lt[::-1], zr[::-1]) * scale

    keep = h < h_cut                                # planes below the cut
    k_above = np.argmin(np.where(h >= h_cut, h, np.inf))
    k_below = np.argmax(np.where(keep, h, -np.inf))
    w = (h_cut - h[k_below]) / (h[k_above] - h[k_below])
    lin = lambda a: (1 - w) * a[:, :, k_below] + w * a[:, :, k_above]
    loglin = lambda a: np.exp(lin(np.log(a)))
    top = [loglin(ne), loglin(temp), lin(vx), lin(vy), lin(vz), loglin(rho)]
    z_top = z[0] if top_first else z[-1]
    z_cut = bottom + np.sign(z_top - bottom) * h_cut

    idx = np.flatnonzero(keep)
    arrays = []
    for a, t in zip((ne, temp, vx, vy, vz, rho), top):
        kept = a[:, :, idx]
        arrays.append(np.concatenate([t[:, :, None], kept] if top_first else [kept, t[:, :, None]], axis=2))
    z_new = np.r_[z_cut, z[idx]] if top_first else np.r_[z[idx], z_cut]

    out.mkdir(parents=True, exist_ok=True)
    with (out / snap.name).open('wb') as f:
        for a in arrays:
            a.astype(np.float32).ravel(order='F').tofile(f)
    write_mesh(out / mesh.name, (x, y, z_new))
    with (out / 'crop.txt').open('a') as f:
        f.write(f'{snap} {out / snap.name} {h_cut:.6e} {H:.6e} {ltau:.2f} {len(z_new)}\n')
    print(f'{snap.name}: cut at {h_cut / 1e5:.0f} of {H / 1e5:.0f} km (<log tau_500> = {ltau}), '
          f'{len(z_new)} of {nz} planes kept; box top <log tau_500> was {lt[-1]:.2f}')


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], *map(float, sys.argv[4:5]))
