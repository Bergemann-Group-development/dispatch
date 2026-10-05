"""The crop agreement test (make_nml.py crop_test, ctest.sbatch): cropped against uncropped
3D Fe departures, and the stage-2 SEDs on them.

    python ctest_compare.py [RUN_DIR]        # default runs/ctest/t60g40m20_s00091

Departures: log b of fe_full against dep_crop (fe_crop's departures mapped in height onto
the full mesh, map_departures.py), per Fe level and cell, binned in the plane mean of
log tau_500 of the full mesh.  Above the cut dep_crop holds the cropped run's top value,
so only the bins below it test the crop.

SEDs, per LINES_SED line, as abundance differences corr(a, b) = A(a) - A(b) at the same
EW with the slope of the nearest batch-1a run (delta_maps.py):
    crop     sed_crop  - sed_full        (NLTE; the crop's effect on the spectrum)
    head     sed_full_head - sed_full    (NLTE; 4890f9fde against the production build)
    lte_sed  sed_full_head LTE set - lte (does the HEAD build's LTE set have the NLTE
                                          species' lines?  Same for sed_full's)
"""
import sys
from pathlib import Path

import numpy as np

import tumult
from delta_maps import RUNS, LINES_SED, ew, slopes, nearest_b1, corr

CUT = -5.0


def logb(run, folder, element='Fe'):
    pop = np.asarray(run.read_patch_save(f'pop_{element}', folder=Path(folder) / 'save', zdim=2))
    return np.log10(np.maximum(pop[:, :, :, :, 1, 0], 1e-30))      # (nx, ny, nz, level)


def departures(top):
    full = tumult.read(str(top / 'fe_full'))
    a = logb(full, top / 'fe_full')
    b = logb(full, top / 'dep_crop')
    ltau = np.asarray(full.atmos.ltau, float)
    plane = ltau.mean(axis=(0, 1))
    d = np.abs(a - b)                                                # (nx, ny, nz, level)
    print(f'Fe departures, |dlog b| fe_full - dep_crop ({a.shape[3]} levels, '
          f'{a.shape[0]}x{a.shape[1]} columns)')
    print(' <log tau_500>    planes  median     p99      max   (level, z of max)')
    edges = [-9, -6, CUT, -4.5, -4, -3, -2, -1, 0, 3]
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = np.where((plane >= lo) & (plane < hi))[0]
        if not len(k):
            continue
        dk = d[:, :, k, :]
        i = np.unravel_index(np.argmax(dk), dk.shape)
        print(f'  [{lo:5.1f},{hi:5.1f})  {len(k):6d}  {np.median(dk):.2e}  '
              f'{np.percentile(dk, 99):.2e}  {dk.max():.2e}   ({i[3] + 1}, {k[i[2]]})')


def seds(top):
    names = sorted(p.name for p in (RUNS / 'b1/hrs').iterdir())
    s = slopes(f'b1/hrs/{nearest_b1(6000, 4.0, -2.0, names)}', LINES_SED)
    runs = {k: tumult.read(str(top / k)) for k in ('sed_full', 'sed_crop', 'sed_full_head', 'lte')}
    print('\nSEDs: A(a) - A(b) at the same EW, dex; EWs in mA')
    print(f'{"line":>12} {"EW full":>8} {"crop":>8} {"head":>8} {"LTE lte":>8} '
          f'{"LTEset full":>11} {"LTEset head":>11}')
    for el, lam0, hw in LINES_SED:
        full = ew(runs['sed_full'], lam0, hw)
        crop = ew(runs['sed_crop'], lam0, hw)
        head = ew(runs['sed_full_head'], lam0, hw)
        lte = ew(runs['lte'], lam0, hw, LTE=True)
        lte_full = ew(runs['sed_full'], lam0, hw, LTE=True)
        lte_head = ew(runs['sed_full_head'], lam0, hw, LTE=True)
        sl = s[lam0]
        print(f'{el:>4} {lam0:7.1f} {full:8.1f} {corr(crop, full, sl):+8.4f} '
              f'{corr(head, full, sl):+8.4f} {lte:8.1f} {lte_full:11.1f} {lte_head:11.1f}')


if __name__ == '__main__':
    top = RUNS / (sys.argv[1] if len(sys.argv) > 1 else 'ctest/t60g40m20_s00091')
    departures(top)
    seds(top)
