from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import numpy as np
import pandas as pd
from pathlib import Path
from tumult.auxfuncs import update_dict
import tumult as tm


def _run_single(nml, da, name, ref, top, launch_kwargs):
    """Top-level function so it's picklable for ProcessPoolExecutor."""
    nml = update_dict(nml, {'atom_params': {'abundance': da}})
    if ref == 'NLTE' and da != 0.0: 
      nml = update_dict(nml, {'m3d_params': {'maxiter': 0}})
    nml_name = Path(top) / f'{name}_{da*100:+04.0f}'
    out, ok = tm.launch(nml_name, nml, **launch_kwargs)
    if not ok:
        print(f'Run {nml_name} failed!')
        return da, None
    return da, tm.read(out)


def _run_batch(da_values, nml, name, ref, top, launch_kwargs, max_threads):
    """Run multiple da values in parallel, splitting threads across workers."""
    n_runs = len(da_values)
    if n_runs == 0:
        return {}

    # Each worker gets a fair share; at least 1 thread each
    threads_per_run = max(1, max_threads // n_runs)
    n_workers = min(n_runs, max_threads)  # no more workers than threads

    worker_kwargs = dict(launch_kwargs, threads=threads_per_run)

    results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_run_single, nml, da, name, ref, top, worker_kwargs): da
            for da in da_values
        }
        for future in as_completed(futures):
            da, result = future.result()
            results[da] = result

    return results


def get_NLTE_corrections(top, stars, atom_file, lines, atmos_files, absmet_files=None, settings=None,
                         tol=0.01, limit=1.0, ref='LTE', integrate=True, verbose=False,
                         max_threads=None, launch_kwargs=None, crop_kwargs=None):

    launch_kwargs = launch_kwargs if launch_kwargs is not None else {}
    crop_kwargs = crop_kwargs if crop_kwargs is not None else {}

    if max_threads is None:
        max_threads = os.cpu_count() or 1

    # Remove 'threads' from launch_kwargs — _run_batch manages it
    launch_kwargs = {k: v for k, v in launch_kwargs.items() if k != 'threads'}

    df = pd.DataFrame()
    all_runs = {}
    atom_file = Path(atom_file)
    sign = 1 if ref == 'LTE' else -1
    max_iter = 10

    for i, star in enumerate(stars):
        teff, logg, feh = star[:3]
        vmic = star[3] if len(star) == 4 else 0.0
        df.loc[i, ['Teff', 'logg', '[Fe/H]', 'vmic'][:len(star)]] = star

        name = atom_file.stem if atom_file.stem != 'atom' else atom_file.suffix[1:]
        name += f'_t{teff/100:.0f}g{logg*10:02.0f}m{feh*10:+03.0f}'
        if vmic > 0:
            name += f'v{vmic:.0f}'

        if star in atmos_files:
            atmos_file = atmos_files[star]
        else:
            print(f'No atmosphere file for star {star}!')
            continue

        if feh in absmet_files:
            absmet_file = absmet_files[feh]
        else:
            print(f'No absmet file for [Fe/H]={feh}!')
            absmet_file = ''

        nml = {
            'timer_params': {'sec_per_report': 1e8, 'verbose': 2},
            'atmos_params': {'atmos_file': atmos_file, 'feh': feh, 'vmic': vmic},
            'atom_params': {
                'atom_file': atom_file,
                'relative_abundance': True,
                'use_atom_abnd': False,
                'convlim': 0.01,
                'convmax': 0.001,
            },
            'm3d_params': {'maxiter': 99, 'ilambd': 3},
            'composition_params': {'absmet_file': absmet_file},
        }
        if settings is not None:
            nml = update_dict(nml, settings)

        def build_grid():
            if any(runs[d] is None for d in runs):
                raise RuntimeError("Some runs failed, cannot build grid.")
            das = sorted(runs.keys())
            return tm.grid([runs[d] for d in das], da=das), das

        def fit_all_lines(grid):
            guesses = {}
            for l0 in lines:
                line = grid.spectra[0].get_line(l0, verbose=False)
                kw = dict(crop_kwargs, xmin=line.lb, xmax=line.lr)
                pfit, _ = grid.fit(refs[l0], fit=['da'], LTE=(ref != 'LTE'), integrate=integrate, **kw)
                guesses[l0] = np.clip(sign * pfit['da'], -limit, limit)
            return guesses

        def find_interval(da_val, knots):
            idx = np.searchsorted(knots, da_val, side='right') - 1
            return max(0, min(idx, len(knots) - 2))

        # Step 1: Build initial coarse grid — all 3 runs in parallel
        da_grid = sorted([-limit, 0.0, limit])
        runs = _run_batch(da_grid, nml, name, ref, top, launch_kwargs, max_threads)
        grid, das = build_grid()

        refs = {}
        for l0 in lines:
          line = grid.spectra[0].get_line(l0, verbose=i==0)
          kw = dict(crop_kwargs, xmin=line.lb, xmax=line.lr)
          refs[l0] = grid.crop(da=0.0, LTE=(ref == 'LTE'), **kw)
        
        guesses = fit_all_lines(grid)

        for iteration in range(max_iter):
            prev_guesses = dict(guesses)
            knots = np.array(sorted(runs.keys()))

            intervals = defaultdict(list)
            for l0, g in guesses.items():
                ii = find_interval(g, knots)
                intervals[(ii, knots[ii], knots[ii + 1])].append((l0, g))

            new_points = []
            for (ii, lo, hi), line_guesses in intervals.items():
                if hi - lo < tol:
                    continue
                guess_values = [g for _, g in line_guesses]
                split_at = np.median(guess_values)

                # Snap to nearest multiple of tol
                split_at = round(split_at / tol) * tol
            
                if (split_at - lo) < tol or (hi - split_at) < tol:
                    split_at = (lo + hi) / 2.0
                if (split_at - lo) < tol or (hi - split_at) < tol:
                    continue

                # Avoid exact duplicates of existing grid points
                if split_at not in runs:
                    new_points.append(split_at)

            if not new_points:
                break

            # Run all new points in parallel
            new_runs = _run_batch(new_points, nml, name, ref, top, launch_kwargs, max_threads)
            runs.update(new_runs)

            grid, das = build_grid()
            guesses = fit_all_lines(grid)

            max_shift = max(abs(guesses[l0] - prev_guesses[l0]) for l0 in lines)
            
            if verbose:
              print(f"Iter {iteration}: added {len(new_points)} points, "
                    f"grid size={len(runs)}, max shift={max_shift:.2e}")

            if max_shift < tol:
                break

        for l0 in lines:
            df.loc[i, l0] = guesses[l0]

        all_runs[star] = runs

        if i==0:
          print(df.to_string(float_format="{:.2f}".format))
        else:
          print('\n'.join(df.loc[i:i+1].to_string(float_format="{:.2f}".format).split('\n')[1:]))

    return df, all_runs