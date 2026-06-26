from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import hashlib
import itertools
import json
import math
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


def _output_folder(nml_name, nml):
    """Where launch() will write this job's run (mirrors launch.py's path logic)."""
    datadir = nml.get('io_params', {}).get('datadir', 'data')
    if isinstance(datadir, str):
        datadir = datadir.strip()
    return Path(datadir) / nml_name


def _run_batch(jobs, launch_kwargs, max_threads):
    """Run jobs, returning {key: run_or_None}.

    Cache hits (output already on disk) are just file reads, so do them inline; only jobs
    that must actually run dispatch.x go to the process pool — spawning workers (each
    re-imports the scientific stack) costs seconds and is pure waste for cache reads.
    """
    if not jobs:
        return {}

    results = {}
    todo = []
    for key, nml_name, nml in jobs:
        if _output_folder(nml_name, nml).exists() and not launch_kwargs.get('overwrite', False):
            results[key] = _run_single(key, nml_name, nml,
                                       dict(launch_kwargs, threads=1))[1]
        else:
            todo.append((key, nml_name, nml))

    if todo:
        threads_per_run = max(1, max_threads // len(todo))
        n_workers = min(len(todo), max_threads)
        worker_kwargs = dict(launch_kwargs, threads=threads_per_run)
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_run_single, key, nml_name, nml, worker_kwargs): key
                for key, nml_name, nml in todo
            }
            for future in as_completed(futures):
                key, result = future.result()
                results[key] = result

    return results


def _vmic_levels(vmic, step):
    """Lattice points bracketing ``vmic`` on a ``step``-spaced grid, with weights.

    Returns ``[(v, w), ...]`` summing to 1. ``vmic`` is synthesised by overriding the
    atmosphere's microturbulence (a uniform Doppler term, cheap to vary on a fixed
    atmosphere), so we snap to a regular lattice — reused across nearby queries via the
    run cache — and interpolate linearly between the two enclosing lattice points.
    """
    vmic = max(0.0, float(vmic))
    lo = math.floor(vmic / step) * step
    hi = math.ceil(vmic / step) * step
    if math.isclose(lo, hi):
        return [(lo, 1.0)]
    alpha = (hi - vmic) / (hi - lo)
    return [(lo, alpha), (hi, 1.0 - alpha)]


def _is_nlte(settings):
    """True when ``settings`` describe an NLTE solve (an atom plus iterations)."""
    return ('atom_params' in (settings or {})
            and (settings.get('m3d_params', {}) or {}).get('maxiter', 0) > 0)


def build_stellar_grid(top, star, atmos_files,
                       settings=None, absmet_files=None,
                       method='Delaunay', max_threads=None, launch_kwargs=None,
                       vmic_step=0.5, departures_top=None):
    """
    Build a tm.grid for the target stellar point by finding the enclosing
    model atmospheres, launching MULTI3D runs for them in parallel, and
    loading the results into a SpectrumGrid indexed by stellar parameters.

    vmic is *decoupled* from atmosphere selection: it barely affects the atmospheric
    structure but strongly affects line strengths, and MULTI3D applies it as a uniform
    namelist override. So atmospheres are bracketed only in (teff, logg, feh); for each
    selected node we take the single atmosphere nearest the query vmic and synthesise it
    at the two ``vmic_step``-spaced lattice points bracketing the query, then interpolate
    the spectra in vmic. This gives fine vmic sampling from a coarse atmosphere grid and
    lets vmic exceed the grid's range.

    Parameters
    ----------
    top : str or Path
        Output directory for the runs.
    star : tuple
        Target stellar parameters as (teff, logg, feh) or (teff, logg, feh, vmic).
    atmos_files : dict
        {(teff, logg, feh) or (teff, logg, feh, vmic): path} mapping. When keys carry a
        vmic, it is used only to pick the representative atmosphere per (teff, logg, feh).
    settings : dict, optional
        NameList settings (atom_params, m3d_params, …) merged into each
        per-atmosphere nml. atmos_params and composition_params are filled
        automatically per atmosphere.
    absmet_files : dict, optional
        {feh: path} mapping for absmet files.
    method : {'Delaunay', 'NaturalNeighbor'}
        How to select the enclosing atmosphere nodes (in teff/logg/feh).
    max_threads, launch_kwargs : optional
        Forwarded to _run_batch / tm.launch.
    vmic_step : float
        Microturbulence lattice spacing in km/s (default 0.5).
    departures_top : str or Path, optional
        When given *and* ``settings`` describe an NLTE solve, run NLTE in two steps: first
        solve departure coefficients (no line list, no output spectrum) for each atmosphere
        under ``departures_top``, saving them via the 'save' atmosphere format; then read
        them back in the line-synthesis step. ``departures_top`` is kept independent of the
        window/line list by the caller so departures are reused across them. When omitted,
        NLTE runs as a single combined solve+synthesis (legacy behaviour).

    Returns
    -------
    tm.grid with a .target attribute {'teff', 'logg', 'feh', 'vmic'}.
    """
    teff, logg, feh = star[:3]
    vmic = star[3] if len(star) == 4 else 0.0

    launch_kwargs = {k: v for k, v in (launch_kwargs or {}).items() if k != 'threads'}
    if max_threads is None:
        max_threads = os.cpu_count() or 1

    # Collapse the atmosphere grid to unique (teff, logg, feh): vmic is no longer an
    # interpolation axis, so for each spatial node we keep one representative atmosphere —
    # the one whose grid vmic is closest to the query (its structure best matches; vmic is
    # overridden anyway). atmos keys may be 3- or 4-tuples.
    rep = {}  # (teff, logg, feh) -> (grid_vmic, path)
    for key, path in atmos_files.items():
        spatial = tuple(float(x) for x in key[:3])
        kvmic = float(key[3]) if len(key) == 4 else 0.0
        cur = rep.get(spatial)
        if cur is None or abs(kvmic - vmic) < abs(cur[0] - vmic):
            rep[spatial] = (kvmic, path)

    spatial_keys = list(rep.keys())
    points = np.array(spatial_keys, dtype=float)
    query = np.array([teff, logg, feh])

    # select interpolation nodes (only in the dimensions where the query is off-grid)
    from tumult.spectrum_grid import get_interp_weights
    weights, indices = get_interp_weights(points, query, method)

    levels = _vmic_levels(vmic, vmic_step)
    base = copy.deepcopy(settings or {})
    nlte = departures_top is not None and _is_nlte(base)

    def run_name(teff_i, logg_i, feh_i, v):
        return f't{teff_i:04.0f}g{logg_i*10:02.0f}m{feh_i*10:+03.0f}v{round(v*10):02.0f}'

    # Cartesian product of selected spatial nodes x vmic lattice points. A node's combined
    # weight is the spatial node weight times the vmic interpolation weight; the same
    # atmosphere may run at two vmic levels (distinct keys).
    nodes = []  # (jkey, teff_i, logg_i, feh_i, v, weight, atmos_path, absmet_file)
    for w_spatial, idx in zip(weights, indices):
        spatial = spatial_keys[idx]
        teff_i, logg_i, feh_i = spatial
        atmos_path = rep[spatial][1]
        absmet_file = absmet_files.get(feh_i, '') if absmet_files else ''
        for v, w_vmic in levels:
            nodes.append(((spatial, v), teff_i, logg_i, feh_i, v,
                          float(w_spatial) * float(w_vmic), atmos_path, absmet_file))

    # NLTE step 1: solve departure coefficients per atmosphere (no line list, no output
    # spectrum), saving the atmosphere + departures via the 'save' format. These runs live
    # under `departures_top` (keyed independently of window/line list), so they are reused
    # across windows and line lists. The synthesis step below reads them back.
    dep_dirs = {}
    if nlte:
        dep_jobs = []
        for jkey, teff_i, logg_i, feh_i, v, _w, atmos_path, absmet_file in nodes:
            dep_nml = copy.deepcopy(base)
            dep_nml.pop('linelist_params', None)
            dep_nml.pop('spectrum_params', None)
            dep_nml = update_dict(dep_nml, {
                'atmos_params': {'atmos_file': atmos_path, 'feh': feh_i, 'vmic': v},
                'composition_params': {'absmet_file': absmet_file},
                # only the departures are wanted here; skip the (expensive) formal solution
                # for an output spectrum (long_scheme='none' => Spectrum.has_spectrum False).
                'm3d_params': {'maxiter': 99, 'ilambd': 3, 'long_scheme': 'none'},
            })
            dep_jobs.append((jkey, Path(departures_top) / run_name(teff_i, logg_i, feh_i, v),
                             dep_nml))
        dep_runs = _run_batch(dep_jobs, launch_kwargs, max_threads)
        dep_failed = [n[0] for n in nodes if dep_runs.get(n[0]) is None]
        if dep_failed:
            raise RuntimeError(f"Departure runs failed for atmospheres: {dep_failed}")
        dep_dirs = {n[0]: str(dep_runs[n[0]].folder) for n in nodes}

    # Synthesis: LTE reads the marcs atmosphere directly; NLTE reads the saved atmosphere +
    # precomputed departures ('save' format, maxiter=0 so no re-solve).
    jobs = []
    job_keys, job_weights, job_params = [], [], []
    for jkey, teff_i, logg_i, feh_i, v, weight, atmos_path, absmet_file in nodes:
        if nlte:
            atmos_override = {'atmos_format': 'save', 'atmos_file': dep_dirs[jkey],
                              'read_departures': True, 'feh': feh_i, 'vmic': v}
            extra = {'m3d_params': {'maxiter': 0}}
        else:
            atmos_override = {'atmos_file': atmos_path, 'feh': feh_i, 'vmic': v}
            extra = {}
        nml_i = update_dict(copy.deepcopy(base), {
            'atmos_params': atmos_override,
            'composition_params': {'absmet_file': absmet_file},
            **extra,
        })
        jobs.append((jkey, Path(top) / run_name(teff_i, logg_i, feh_i, v), nml_i))
        job_keys.append(jkey)
        job_weights.append(weight)
        job_params.append((teff_i, logg_i, feh_i, v))

    runs = _run_batch(jobs, launch_kwargs, max_threads)

    failed = [k for k in job_keys if runs.get(k) is None]
    if failed:
        raise RuntimeError(f"Runs failed for atmospheres: {failed}")

    param_dict = {
        'teff': [p[0] for p in job_params],
        'logg': [p[1] for p in job_params],
        'feh':  [p[2] for p in job_params],
        'vmic': [p[3] for p in job_params],
    }

    g = tm.grid([runs[k] for k in job_keys], **param_dict)
    g.p0 = {'teff': teff, 'logg': logg, 'feh': feh, 'vmic': vmic}
    # the interpolation weights for the selected nodes (aligned with g.spectra); the
    # caller can combine directly with these instead of re-deriving them on the subset.
    g.interp_weights = np.asarray(job_weights, dtype=float)
    return g


def compute_spectrum(top, star, atmos_files, *, settings=None, absmet_files=None,
                     method='Delaunay', max_threads=None, launch_kwargs=None,
                     xmin=None, xmax=None, pad=None, vmic_step=0.5):
    """Compute a single spectrum at an arbitrary stellar point (Teff, logg, FeH, vmic).

    Thin wrapper over :func:`build_stellar_grid`: it synthesises the enclosing model
    atmospheres in parallel and interpolates the result at ``star``. LTE vs NLTE is
    encoded entirely in ``settings`` (presence of ``atom_params`` / ``m3d_params.maxiter``);
    this helper only fills in the synthesis wavelength window from ``xmin``/``xmax`` so the
    Fortran code computes just the visible range.

    Parameters
    ----------
    top, star, atmos_files, settings, absmet_files, method, max_threads, launch_kwargs
        Forwarded to :func:`build_stellar_grid` (``settings`` is the fully assembled nml
        dict, keyed by namelist section).
    xmin, xmax : float, optional
        Wavelength window in Angstrom. Sets ``spectrum_params.aa_blue/aa_red`` (unless a
        ``lam_file`` is configured) and crops the interpolated spectrum to this range.
    pad : float, optional
        Extra Angstrom synthesised on each side of [xmin, xmax] so the cropped window has
        clean edges for convolution. Defaults to ``max(0.5, 0.05 * (xmax - xmin))``.

    Returns
    -------
    Spectrum
        The continuum-normalised, interpolated spectrum at ``star``.
    """
    settings = copy.deepcopy(settings or {})
    top = Path(top)
    label = 'full'

    if xmin is not None and xmax is not None:
        lo, hi = (xmin, xmax) if xmin <= xmax else (xmax, xmin)
        sp = settings.setdefault('spectrum_params', {})
        if not sp.get('lam_file'):
            lo_i, hi_i = int(np.floor(lo)), int(np.ceil(hi))
            if pad is None:
                pad = min(max(0.5, 0.05 * (hi_i - lo_i)), 100)
            sp.setdefault('aa_blue', float(lo_i - pad))
            sp.setdefault('aa_red', float(hi_i + pad))
            label = f'w{lo_i}_{hi_i}'

    # build_stellar_grid names cached runs by stellar params only, so a run computed with
    # different settings (window, R, abundances, mode, line list, ...) would be silently
    # reused. Key the output path on a hash of the settings so any change re-synthesises;
    # identical settings still hit the disk cache.
    digest = hashlib.sha1(
        json.dumps(settings, sort_keys=True, default=str).encode()).hexdigest()[:8]
    syn_top = top / f'{label}_{digest}'

    # NLTE departure runs are independent of the window and line list, so key their cache on
    # only the departure-relevant settings and place them in a sibling `dep_<hash>` bucket —
    # reused across windows/line lists (see build_stellar_grid's two-step path).
    departures_top = None
    if _is_nlte(settings):
        dep_subset = {k: settings[k] for k in
                      ('atom_params', 'composition_params', 'atmos_params', 'm3d_params')
                      if k in settings}
        dep_digest = hashlib.sha1(
            json.dumps(dep_subset, sort_keys=True, default=str).encode()).hexdigest()[:8]
        departures_top = top / f'dep_{dep_digest}'

    grid = build_stellar_grid(str(syn_top), star, atmos_files, settings=settings,
                              absmet_files=absmet_files, method=method,
                              max_threads=max_threads, launch_kwargs=launch_kwargs,
                              vmic_step=vmic_step,
                              departures_top=str(departures_top) if departures_top else None)

    # a run can finish yet emit no spectrum (e.g. no R/daa/nnu/lam_file set, or
    # long_scheme='none'); flag it clearly instead of failing later on a missing attribute.
    no_spec = [s for s in grid.spectra if not getattr(s, 'has_spectrum', False)]
    if no_spec:
        folders = ', '.join(str(getattr(s, 'folder', '?')) for s in no_spec)
        raise RuntimeError(
            f"Synthesis produced no spectrum for: {folders}. Check the spectrum "
            f"settings (R / wavelength range) and the synthesis line list.")

    # combine the selected nodes using the (truncated) selection weights — the same
    # weighted average grid.crop would do, but reusing the full-grid weights directly.
    # Pass LTE explicitly per run so an LTE run yields its (only) LTE spectrum without the
    # "Cannot get NLTE spectrum" warning Run.get_spectrum prints on the NLTE default.
    weights = getattr(grid, 'interp_weights', None)
    if weights is not None:
        spec = tm.Spectrum.combine(
            [s.crop(xmin=xmin, xmax=xmax, LTE=(getattr(s, 'mode', 'LTE') == 'LTE'))
             for s in grid.spectra], weights=weights)
    else:
        lte = all(getattr(s, 'mode', 'LTE') == 'LTE' for s in grid.spectra)
        spec = grid.crop(xmin=xmin, xmax=xmax, method=method, LTE=lte)
    if len(spec) == 0:
        raise RuntimeError(
            f"Interpolated spectrum is empty over [{xmin}, {xmax}] Å — the synthesis "
            f"did not cover the requested window.")
    return spec


def adaptive_fit(target, synth_point_fn, axes, guesses, *, post_fit=(), spacings=None,
                 axis_nodes=None, atmos_bounds=None, max_rounds=6, tol=None, xtol=None,
                 fit_kwargs=None):
    """Fit ``target`` by adaptive grid synthesis.

    Builds a minimal grid (2 points per fitted axis, bracketing the guess), synthesises each
    point via ``synth_point_fn``, and fits every axis at once with :meth:`SpectrumGrid.fit`
    (which also fits the continuous ``post_fit`` params — vshift/vrot/zeta_rt — exactly, with
    no grid axis). It then expands any axis whose optimum sits on a boundary, or subdivides
    around an interior optimum that is still coarsely sampled, and refits — until every axis
    optimum is interior and resolved to ``tol`` (or ``max_rounds`` is reached).

    Parameters
    ----------
    target : Spectrum
        Observed spectrum to match (compared on the synthetic grid via ``crop(newx=...)``).
    synth_point_fn : callable
        ``synth_point_fn(point: dict) -> Spectrum`` synthesising a normalised spectrum at one
        grid point (keys are ``axes``). May return ``None`` for an un-synthesisable point.
        Memoised here on the rounded point, so repeated points cost nothing.
    axes : sequence of str
        Names of the fitted grid dimensions (e.g. 'teff', 'logg', 'feh', 'vmic', 'Mg').
    guesses : dict
        Initial value for every fitted name — grid ``axes`` *and* ``post_fit`` params.
    post_fit : sequence of str, optional
        Continuous post-processing params (vshift/vrot/zeta_rt) fit by ``Spectrum.fit``.
    spacings : dict, optional
        Step for each *continuous* axis (vmic, abundances) — used to bracket the guess and
        expand outward, with sub-step refinement around an interior optimum. Default 1.0.
    axis_nodes : dict, optional
        ``{axis: sorted values}`` for *discrete* axes whose sample points are fixed (the
        atmosphere-grid nodes for teff/logg/feh). These bracket the guess and expand to
        adjacent nodes but are never sub-sampled (finer than a node is just interpolation
        between the same two syntheses).
    atmos_bounds : dict, optional
        ``{axis: (lo, hi)}`` hard clamp for continuous-axis expansion (``None`` = unbounded),
        e.g. ``(0, None)`` for vmic. Discrete axes are bounded by their node list.
    max_rounds : int, optional
    tol : dict, optional
        Per-axis convergence tolerance for grid expansion/subdivision (default:
        ``spacings[axis] / 4``). This governs how finely the *grid* is sampled, not the
        optimiser precision (see ``xtol``).
    xtol : dict, optional
        Per-parameter absolute optimiser tolerance, forwarded to ``Spectrum.fit`` (e.g.
        ``{"Mg": 0.01, "vshift": 0.1}``). Also the precision the returned best-fit values are
        rounded to. Velocities default to 0.1; everything else to 0.01.
    fit_kwargs : dict, optional
        Extra kwargs forwarded to ``grid.fit`` (e.g. ``method``, ``integrate``, ``sigma``).

    Returns
    -------
    (best : dict, fitted : Spectrum, grid : SpectrumGrid)
    """
    axes = list(axes)
    post_fit = list(post_fit)
    spacings = dict(spacings or {})
    axis_nodes = {k: sorted(float(x) for x in v) for k, v in (axis_nodes or {}).items()}
    atmos_bounds = dict(atmos_bounds or {})
    fit_kwargs = dict(fit_kwargs or {})
    if xtol:
        fit_kwargs.setdefault("xtol", dict(xtol))   # optimiser precision per parameter
    tol = dict(tol or {})

    def _clamp(axis, v):
        lo, hi = atmos_bounds.get(axis, (None, None))
        if lo is not None:
            v = max(v, lo)
        if hi is not None:
            v = min(v, hi)
        return v

    def _bracket_nodes(nodes, g):
        """The two nodes bracketing ``g`` (or the nearest two if ``g`` is off the ends)."""
        i = int(np.searchsorted(nodes, g))
        i = max(1, min(i, len(nodes) - 1))
        return [nodes[i - 1], nodes[i]]

    cache = {}

    def _synth(point):
        key = tuple(sorted((k, round(float(v), 6)) for k, v in point.items()))
        if key not in cache:
            cache[key] = synth_point_fn(dict(point))
        return cache[key]

    # initial samples per axis: discrete axes bracket the guess at their nodes; continuous
    # axes take 2 points bracketing the guess (re-spread if clamping collapses them).
    samples = {}
    for a in axes:
        g = float(guesses[a])
        if a in axis_nodes:
            nodes = axis_nodes[a]
            samples[a] = _bracket_nodes(nodes, g) if len(nodes) >= 2 else list(nodes) or [g]
        else:
            s = float(spacings.get(a, 1.0))
            pts = sorted({_clamp(a, g - s / 2), _clamp(a, g + s / 2)})
            if len(pts) < 2:
                base = pts[0] if pts else g
                pts = sorted({base, _clamp(a, base + s), _clamp(a, base - s)})
            samples[a] = pts

    best = fitted = grid = None
    for _round in range(max_rounds):
        # synthesise the Cartesian product of the current samples
        spectra, params = [], {a: [] for a in axes}
        for combo in itertools.product(*(samples[a] for a in axes)):
            point = dict(zip(axes, combo))
            spec = _synth(point)
            if spec is None:
                continue
            spectra.append(spec)
            for a, v in point.items():
                params[a].append(v)
        # need >=2 points to interpolate a grid axis; a post-params-only fit (no axes) is
        # a single synthesis that Spectrum.fit transforms in place.
        if len(spectra) < (2 if axes else 1):
            raise RuntimeError("adaptive_fit: too few synthesised grid points to fit.")

        grid = tm.grid(spectra, **params)
        p0 = {k: float(guesses[k]) for k in axes + post_fit if k in guesses}
        best, fitted = grid.fit(target, fit=axes + post_fit, p0=p0, **fit_kwargs)

        # per-axis: expand at a boundary optimum; for continuous axes also subdivide a coarse
        # interior bracket. Discrete (node) axes only expand to the next node.
        changed = False
        for a in axes:
            pts = samples[a]
            lo, hi = pts[0], pts[-1]
            v = float(best[a])
            if a in axis_nodes:
                nodes = axis_nodes[a]
                t = float(tol.get(a, 0.1 * (hi - lo if hi > lo else 1.0)))
                if v <= lo + t:                       # extend to the previous node
                    i = int(np.searchsorted(nodes, lo))
                    if i > 0 and nodes[i - 1] not in pts:
                        pts.insert(0, nodes[i - 1])
                        changed = True
                elif v >= hi - t:                     # extend to the next node
                    i = int(np.searchsorted(nodes, hi))
                    if i < len(nodes) - 1 and nodes[i + 1] not in pts:
                        pts.append(nodes[i + 1])
                        changed = True
            else:
                s = float(spacings.get(a, 1.0))
                t = float(tol.get(a, s / 4))
                if v <= lo + t:
                    nv = _clamp(a, lo - s)
                    if nv < lo - 1e-9:
                        pts.append(nv)
                        changed = True
                elif v >= hi - t:
                    nv = _clamp(a, hi + s)
                    if nv > hi + 1e-9:
                        pts.append(nv)
                        changed = True
                else:
                    left = max(p for p in pts if p <= v)
                    right = min(p for p in pts if p >= v)
                    if (right - left) > 2 * t and all(abs(v - p) > 1e-9 for p in pts):
                        pts.append(v)
                        changed = True
            pts.sort()

        if not changed:
            break

    return best, fitted, grid


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
    name += f'_t{teff:04.0f}g{logg*10:02.0f}m{feh*10:+03.0f}'
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
