from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import os
import numpy as np
import pandas as pd
from pathlib import Path
from tumult.auxfuncs import update_dict
import tumult as tm


def _run_single(key, nml_name, nml, launch_kwargs):
    """Picklable worker: launch one job and return (key, run_or_None)."""
    out, ok = tm.launch(nml_name, nml, **launch_kwargs)
    if not ok:
        if out.exists():
            return key, tm.read(out)  # run already completed, just read it
        print(f'Run {nml_name} failed!')
        return key, None
    return key, tm.read(out)


def _make_da_jobs(da_values, nml, name, ref, top):
    """Build a jobs list for _run_batch from abundance offsets."""
    jobs = []
    for da in da_values:
        nml_i = update_dict(copy.deepcopy(nml), {'atom_params': {'abundance': da}})
        if ref == 'NLTE' and da != 0.0:
            nml_i = update_dict(nml_i, {'m3d_params': {'maxiter': 0}})
        nml_name = Path(top) / f'{name}_{da*100:+04.0f}'
        jobs.append((da, nml_name, nml_i))
    return jobs


def _run_batch(jobs, launch_kwargs, max_threads):
    """Run jobs in parallel, splitting threads across workers.

    jobs: list of (key, nml_name, nml) triples.
    Returns {key: run_or_None}.
    """
    n_runs = len(jobs)
    if n_runs == 0:
        return {}

    threads_per_run = max(1, max_threads // n_runs)
    n_workers = min(n_runs, max_threads)

    worker_kwargs = dict(launch_kwargs, threads=threads_per_run)

    results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_run_single, key, nml_name, nml, worker_kwargs): key
            for key, nml_name, nml in jobs
        }
        for future in as_completed(futures):
            key, result = future.result()
            results[key] = result

    return results


def build_stellar_grid(top, star, atmos_files,
                       settings=None, absmet_files=None,
                       method='Delaunay', max_threads=None, launch_kwargs=None):
    """
    Build a tm.grid for the target stellar point by finding the enclosing
    model atmospheres, launching MULTI3D runs for them in parallel, and
    loading the results into a SpectrumGrid indexed by stellar parameters.

    Parameters
    ----------
    top : str or Path
        Output directory for the runs.
    star : tuple
        Target stellar parameters as (teff, logg, feh) or (teff, logg, feh, vmic).
    atmos_files : dict
        {(teff, logg, feh) or (teff, logg, feh, vmic): path} mapping.
    settings : dict, optional
        NameList settings (atom_params, m3d_params, …) merged into each
        per-atmosphere nml. atmos_params and composition_params are filled
        automatically per atmosphere.
    absmet_files : dict, optional
        {feh: path} mapping for absmet files.
    method : {'Delaunay', 'NaturalNeighbor'}
        How to select the enclosing atmosphere nodes.
    max_threads, launch_kwargs : optional
        Forwarded to _run_batch / tm.launch.

    Returns
    -------
    tm.grid with a .target attribute {'teff', 'logg', 'feh', 'vmic'}.
    """
    from tumult.spectrum_grid import normalize_points, get_delaunay_weights, get_nn_weights

    teff, logg, feh = star[:3]
    vmic = star[3] if len(star) == 4 else 0.0

    launch_kwargs = {k: v for k, v in (launch_kwargs or {}).items() if k != 'threads'}
    if max_threads is None:
        max_threads = os.cpu_count() or 1

    stars = list(atmos_files.keys())
    vmic_varies = len(stars[0]) == 4

    points = np.array(stars, dtype=float)
    if not vmic_varies:
        points = points[:, :3]
    query = np.array([teff, logg, feh, vmic] if vmic_varies else [teff, logg, feh])

    points_norm, p_min, p_ptp = normalize_points(points)
    query_norm = (query - p_min) / p_ptp

    if method == 'Delaunay':
        _, indices = get_delaunay_weights(points_norm, query_norm)
        if indices is None:
            raise ValueError(f"Target {query} is outside the convex hull of available atmospheres.")
    elif method == 'NaturalNeighbor':
        _, indices = get_nn_weights(points_norm, query_norm)
        if indices[0] < 0:
            raise ValueError(f"Target {query} is outside the convex hull of available atmospheres.")
    else:
        raise ValueError(f"Unknown method: {method!r}. Use 'Delaunay' or 'NaturalNeighbor'.")

    relevant_stars = [stars[i] for i in indices]

    jobs = []
    for star in relevant_stars:
        teff_i, logg_i, feh_i = star[:3]
        vmic_i = float(star[3]) if vmic_varies else vmic
        absmet_file = absmet_files.get(feh_i, '') if absmet_files else ''

        nml_i = update_dict(copy.deepcopy(settings or {}), {
            'atmos_params': {'atmos_file': atmos_files[star], 'feh': feh_i, 'vmic': vmic_i},
            'composition_params': {'absmet_file': absmet_file},
        })

        run_name = f't{teff_i/100:.0f}g{logg_i*10:02.0f}m{feh_i*10:+03.0f}'
        if vmic_i > 0:
            run_name += f'v{vmic_i:.0f}'

        jobs.append((star, Path(top) / run_name, nml_i))

    runs = _run_batch(jobs, launch_kwargs, max_threads)

    failed = [s for s in relevant_stars if runs.get(s) is None]
    if failed:
        raise RuntimeError(f"Runs failed for atmospheres: {failed}")

    param_dict = {
        'teff': [s[0] for s in relevant_stars],
        'logg': [s[1] for s in relevant_stars],
        'feh':  [s[2] for s in relevant_stars],
    }
    if vmic_varies:
        param_dict['vmic'] = [s[3] for s in relevant_stars]

    g = tm.grid([runs[s] for s in relevant_stars], **param_dict)
    g.p0 = {'teff': teff, 'logg': logg, 'feh': feh}
    if vmic_varies:
        g.p0['vmic'] = vmic
    return g


def get_NLTE_correction(top, star, atom_file, atmos_file, lines, absmet_files=None, settings=None,
                        tol=0.01, limit=1.0, ref='LTE', integrate=True, verbose=False,
                        max_threads=None, launch_kwargs=None, crop_kwargs=None):
    """
    Compute NLTE corrections for a single star.

    Returns
    -------
    guesses : dict  {line_center: correction}
    runs    : dict  {da: Run}
    """
    launch_kwargs = {k: v for k, v in (launch_kwargs or {}).items() if k != 'threads'}
    crop_kwargs = crop_kwargs or {}

    if max_threads is None:
        max_threads = os.cpu_count() or 1

    atom_file = Path(atom_file)
    if not atom_file.exists():
        raise FileNotFoundError(f"Atom file {atom_file} not found!")

    sign = 1 if ref == 'LTE' else -1
    max_iter = 10

    teff, logg, feh = star[:3]
    vmic = star[3] if len(star) == 4 else 0.0

    name = atom_file.stem if atom_file.stem != 'atom' else atom_file.suffix[1:]
    name += f'_t{teff/100:.0f}g{logg*10:02.0f}m{feh*10:+03.0f}'
    if vmic > 0:
        name += f'v{vmic:.0f}'

    if absmet_files is not None and feh in absmet_files:
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
    runs = _run_batch(_make_da_jobs(da_grid, nml, name, ref, top), launch_kwargs, max_threads)
    grid, _ = build_grid()

    refs = {}
    for l0 in lines:
        line = grid.spectra[0].get_line(l0, verbose=verbose)
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
        new_runs = _run_batch(_make_da_jobs(new_points, nml, name, ref, top), launch_kwargs, max_threads)
        runs.update(new_runs)

        grid, _ = build_grid()
        guesses = fit_all_lines(grid)

        max_shift = max(abs(guesses[l0] - prev_guesses[l0]) for l0 in lines)

        if verbose:
            print(f"Iter {iteration}: added {len(new_points)} points, "
                  f"grid size={len(runs)}, max shift={max_shift:.2e}")

        if max_shift < tol:
            break

    return guesses, runs


def get_NLTE_corrections(top, stars, atom_file, atmos_files, lines, **kwargs):
    """ Compute NLTE corrections for multiple stars and return a DataFrame.

    Returns
    -------
    df : DataFrame with columns ['Teff', 'logg', '[Fe/H]', 'vmic', line1, line2, ...]
    runs : dict {star: {da: Run}} with the runs for each star.
    """
    atom_file = Path(atom_file)
    if not atom_file.exists():
        raise FileNotFoundError(f"Atom file {atom_file} not found!")

    df = pd.DataFrame()
    runs = {}

    for i, star in enumerate(stars):
        df.loc[i, ['Teff', 'logg', '[Fe/H]', 'vmic'][:len(star)]] = star

        if star not in atmos_files:
            print(f'No atmosphere file for star {star}!')
            continue
        else:
            atmos_file = atmos_files[star]

        corr, runs[star] = get_NLTE_correction(top, star, atom_file, atmos_file, lines, **kwargs)

        for l0 in lines:
            df.loc[i, l0] = corr[l0]

        if i == 0:
            print(df.to_string(float_format="{:.2f}".format))
        else:
            print('\n'.join(df.loc[i:i+1].to_string(float_format="{:.2f}".format).split('\n')[1:]))

    return df, runs