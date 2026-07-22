"""Raw dispatch patch-save reading (binary memmaps + their meta files).

Extracted from the run class in the Phase 6 decomposition: this module knows the
on-disk layout — ``<kind>_meta.txt`` describing per-patch ``.bin`` memmaps keyed by
``(iz, inu)`` — and how to stitch the patches back together. It knows nothing about
the physics of what the arrays mean (see :mod:`tumult.flux` for the intensity layout).
"""
from __future__ import annotations

import os
import pathlib
import textwrap

import numpy as np

from .memmap_list import MemmapList


def load_patch_files(folder, kind):
    """Memmap every patch file of ``kind`` under ``folder`` (a run's ``save/`` dir).

    Returns ``(dct, dims)``: patches keyed by their ``(iz, inu)`` id tuple, plus the
    dimensions of the last-listed patch (the historical default for the concat axis).
    """
    resolved_folder = str(pathlib.Path(folder).resolve()) + os.sep
    meta_file = resolved_folder + f'{kind.replace(" ", "_")}_meta.txt'
    dct = {}
    dims = ()
    try:
        with open(meta_file) as meta:
            for line in meta:
                if line.strip().startswith('*'):
                    continue
                info = line.split()
                file = info[0]
                dtype = f'<f{int(info[1]) // 8}'
                dims = tuple(map(int, info[2:]))
                key = file.split('.')[0].split('_')[-1]
                key = tuple(map(int, textwrap.wrap(key, 3)))
                dct[key] = np.memmap(resolved_folder + file, dtype=dtype, mode='r',
                                     shape=dims, order='F')
    except FileNotFoundError:
        raise FileNotFoundError(f'{meta_file} not found!')
    return dct, dims


def concat_patches(dct, n_nu, fdim, zdim, lazy):
    """Stitch the ``(iz, inu)``-keyed patches into full arrays.

    ``n_nu`` is the run's frequency-block size (``run_params/n_nu``): when patches are
    frequency-split, the positive ``inu`` ids up to ``n_nu``, above ``n_nu``, and the
    negative ids form separate sets (returned as a list, in that order). ``lazy``
    stitches with :class:`MemmapList` (no copy) instead of ``np.concatenate``.
    """
    all_iz = np.unique([key[0] for key in dct.keys()])
    all_inu = np.unique([key[1] for key in dct.keys()])
    f_indep = all(all_inu == 0) or fdim is None

    concat = MemmapList if lazy else np.concatenate

    def merge_set(inus):
        merged = []
        for iz in all_iz:
            parts = [dct[(iz, inu)] for inu in inus if (iz, inu) in dct]
            if not parts:
                continue
            mem = parts[0] if len(parts) == 1 else concat(parts, axis=fdim)
            merged.append(mem)
        if len(merged) == 0:
            return None
        elif len(merged) == 1:
            return merged[0]
        else:
            return concat(merged, axis=zdim)

    if f_indep:
        return merge_set(all_inu)
    else:
        set1 = [inu for inu in all_inu if inu > 0 and inu <= n_nu]
        set2 = [inu for inu in all_inu if inu > 0 and inu > n_nu]
        set3 = [inu for inu in all_inu if inu < 0][::-1]
        sets = [s for s in [set1, set2, set3] if len(s) > 0]
        return [merge_set(s) for s in sets]


def read_patch_save(folder, kind, n_nu, concat=True, fdim=0, zdim=None, lazy=True):
    """Load ``kind`` from a run's save folder: memmap the patches and (by default)
    stitch them together (``zdim`` defaults to the patches' last axis)."""
    dct, dims = load_patch_files(folder, kind)
    if not concat:
        return dct
    zdim = zdim or len(dims) - 1
    return concat_patches(dct, n_nu, fdim, zdim, lazy=lazy)
