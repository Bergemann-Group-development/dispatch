"""Map 1D departure coefficients onto a 3D snapshot in log tau_500, as a read_departures
source for a 3D spectrum (the <3D>-departure approximation of 3D NLTE).

    python map_departures.py RUN_1D RUN_3D OUT ELEMENT [ELEMENT ...]

RUN_1D is a 1D NLTE run of the elements (e.g. the pilot's av3d_mg on the <3D> average
of RUN_3D's model), RUN_3D any 3D run of the snapshot that saved its atmosphere (the
pilot's s<snap>_lte or s<snap>_mg).  Every cell of the snapshot takes, per level,

    log b(cell) = interp(log tau_500(cell), log tau_500(1D), log b(1D)),

linear in log tau_500 and held constant beyond the 1D depth range.  log tau_500 is the
vertical optical depth of each column, as tumult computes it from the saved 500 nm
opacity (slot chi) of each run.

OUT gets params.nml, atmos_mesh.txt and save/atmos_* of RUN_3D as symlinks, and
save/pop_<ELEMENT>_* written in the layout of a 3D NLTE run: one file per z-slab patch,
(nx, ny, nz_patch, nlevel, 2, 1) float32, slot 2 the departures.  Only the departures of
the unrotated mesh are written: read_departures loads that copy and projects the
ray-aligned copies itself (nlte_patch_t%project_populations).  Slot 1 (the LTE
populations of the solving run) is not read by read_departures and is written as 1.
"""
import sys
from pathlib import Path

import numpy as np

import tumult


def profile(run_1d, element):
    """log tau_500 (ascending) and log b (depth, level) of a 1D run."""
    r = tumult.read(str(run_1d))
    pop = np.asarray(r.read_patch_save(f'pop_{element}', folder=Path(run_1d) / 'save', zdim=2))
    b = pop[0, 0, :, :, 1, 0]                       # (z, level), unrotated copy
    ltau = np.asarray(r.atmos.ltau, float).ravel()
    order = np.argsort(ltau)
    return ltau[order], np.log10(np.maximum(b[order], 1e-30))


def write_patches(save, element, b3d, nz_patches):
    """pop_<element>_<iz>000.bin per z-slab and the meta file listing them."""
    rows = ['* filename                     bytes']
    z0 = 0
    for iz, nzp in enumerate(nz_patches, 1):
        part = b3d[:, :, z0:z0 + nzp, :]
        z0 += nzp
        mem = np.ones(part.shape + (2, 1), np.float32)
        mem[..., 1, 0] = part
        name = f'pop_{element}_{iz:03d}000.bin'
        mem.ravel(order='F').tofile(save / name)
        rows.append(f'{name:<28s}{32:>6d}' + ''.join(f'{n:>10d}' for n in mem.shape))
    (save / f'pop_{element}_meta.txt').write_text('\n'.join(rows) + '\n')


def main(run_1d, run_3d, out, elements):
    run_1d, run_3d, out = Path(run_1d), Path(run_3d), Path(out)
    save = out / 'save'
    save.mkdir(parents=True, exist_ok=True)
    for name in ('params.nml', 'atmos_mesh.txt'):
        if not (out / name).is_symlink():
            (out / name).symlink_to(run_3d / name)
    for f in (run_3d / 'save').glob('atmos*'):
        if not (save / f.name).is_symlink():
            (save / f.name).symlink_to(f)
    r3 = tumult.read(str(run_3d))
    ltau3 = np.asarray(r3.atmos.ltau, float)        # (nx, ny, nz), the unrotated mesh
    nz_patches = np.diff(np.r_[0, r3.atmos.atmos_cuts]).astype(int)
    assert nz_patches.sum() == ltau3.shape[2]
    for el in elements:
        lt1, lb1 = profile(run_1d, el)
        nlev = lb1.shape[1]
        b3d = np.empty(ltau3.shape + (nlev,), np.float32)
        for l in range(nlev):
            b3d[..., l] = 10 ** np.interp(ltau3, lt1, lb1[:, l])
        write_patches(save, el, b3d, nz_patches)
        ev = run_1d / 'save' / f'pop_{el}_ev.txt'
        (save / ev.name).write_text(ev.read_text())
        print(f'{el}: {nlev} levels, 1D log tau {lt1[0]:.2f}..{lt1[-1]:.2f}, 3D '
              f'{ltau3.min():.2f}..{ltau3.max():.2f}; b {b3d.min():.3g}..{b3d.max():.3g} -> {save}')


if __name__ == '__main__':
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
