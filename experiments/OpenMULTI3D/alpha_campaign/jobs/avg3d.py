"""<3D> averages from Stagger snapshots, with density.

    python avg3d.py OUT.atmos FEH RUN [RUN ...]

Each RUN is a dispatch.x run on one snapshot whose saved atmosphere carries T, n_e,
rho, the velocities and the 500 nm opacity (any 3D run of the snapshot; make_nml.py
avg3d writes cheap atmosphere-only runs).  The columns of all snapshots are pooled and
averaged on surfaces of constant log tau_500 with tumult's AtmosClass.to_dispatch1d_atmos:
geometric means of T, n_e and rho, and vmic from the velocity dispersion.  The result
is dispatch.x's 'text' format (log tau_500, T, Pe, rho, vmic), which a run reads with
atmos_format='text'; use_ne=F makes it take rho and recompute n_e from (T, rho).

The MULTI-format Stagger <3D> files have no density, and dispatch.x rebuilds it from the
averaged (T, n_e); in metal-poor atmospheres that compresses the tau scale above
log tau -2 (phase 3).  These averages carry the density instead.
"""
import sys
from pathlib import Path

import numpy as np
import tumult
from tumult.atmos import AtmosClass


def pooled(runs):
    """One AtmosClass holding the columns of every run side by side along x."""
    atm = [tumult.read(r).atmos for r in runs]
    out = AtmosClass()
    out.dim = '3D'
    for key in ('temp', 'ne', 'rho', 'vx', 'vy', 'vz', 'ltau'):
        setattr(out, key, np.concatenate([np.asarray(getattr(a, key), float) for a in atm], axis=0))
    return out


def main(out, feh, runs):
    pooled(runs).to_dispatch1d_atmos(Path(out), FeH=feh)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
