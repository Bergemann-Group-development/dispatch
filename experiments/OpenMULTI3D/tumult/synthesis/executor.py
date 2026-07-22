"""Parallel dispatch of MULTI3D jobs through a ProcessPoolExecutor.

Extracted from the synthesis monolith in the Phase 6 decomposition: this module knows how
to run a batch of (key, nml_name, nml) jobs — inline cache reads, worker pools, thread
budgeting, cancellation — and nothing about which jobs to build or what to do with the
resulting runs.
"""
from __future__ import annotations

import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import tumult as tm
from tumult.auxfuncs import update_dict
from tumult.launch import run_complete

from .cancel import CancelledError, _init_pool_worker


def _run_single(key, nml_name, nml, launch_kwargs):
    """Picklable worker: launch one job and return (key, run_or_None)."""
    out, ok = tm.launch(nml_name, nml, **launch_kwargs)
    # launch returns ok=False when it served an existing cache hit; either way the run is only
    # usable if it actually finished (wrote abund.txt). A crashed run leaves a partial folder —
    # report it as failed (-> a clear "Runs failed" error) instead of crashing the reader.
    if not run_complete(out):
        print(f'Run {nml_name} failed!')
        return key, None
    return key, tm.read(out)


def _make_da_jobs(da_values, nml, name, ref, top):
    """Build a jobs list for _run_batch from abundance offsets."""
    jobs = []
    for da in da_values:
        nml_i = update_dict(copy.deepcopy(nml), {'atom_params': {'abundance': da}})
        if ref == 'NLTE' and da != 0.0:
            nml_i = update_dict(nml_i, {'rt_iteration_params': {'maxiter': 0}})
        nml_name = Path(top) / f'{name}_{da*100:+04.0f}'
        jobs.append((da, nml_name, nml_i))
    return jobs


def _output_folder(nml_name, nml):
    """Where launch() will write this job's run (mirrors launch.py's path logic)."""
    datadir = nml.get('io_params', {}).get('datadir', 'data')
    if isinstance(datadir, str):
        datadir = datadir.strip()
    return Path(datadir) / nml_name


def _run_batch(jobs, launch_kwargs, max_threads, cancel_token=None):
    """Run jobs, returning {key: run_or_None}.

    Cache hits (output already on disk) are just file reads, so do them inline; only jobs
    that must actually run dispatch.x go to the process pool — spawning workers (each
    re-imports the scientific stack) costs seconds and is pure waste for cache reads.

    ``cancel_token`` (optional): the live pool is registered on it so a cancel can kill the
    worker process groups (and their dispatch.x children); a cancel surfaces here as
    :class:`CancelledError`.
    """
    if not jobs:
        return {}
    if cancel_token is not None:
        cancel_token.check()

    results = {}
    todo = []
    for key, nml_name, nml in jobs:
        if (run_complete(_output_folder(nml_name, nml))
                and not launch_kwargs.get('overwrite', False)):
            results[key] = _run_single(key, nml_name, nml,
                                       dict(launch_kwargs, threads=1))[1]
        else:
            todo.append((key, nml_name, nml))

    if todo:
        threads_per_run = launch_kwargs.get('threads', max(1, max_threads // len(todo)))
        n_workers = min(len(todo), max_threads)
        worker_kwargs = dict(launch_kwargs, threads=threads_per_run)
        with ProcessPoolExecutor(max_workers=n_workers,
                                 initializer=_init_pool_worker) as executor:
            try:
                if cancel_token is not None:
                    cancel_token.set_pool(executor)
                    cancel_token.check()   # cancelled during setup — abort before submitting
                futures = {
                    executor.submit(_run_single, key, nml_name, nml, worker_kwargs): key
                    for key, nml_name, nml in todo
                }
                for future in as_completed(futures):
                    try:
                        key, result = future.result()
                    except BrokenProcessPool:
                        # workers were killed: a cancel, or a genuine crash
                        if cancel_token is not None and cancel_token.cancelled:
                            raise CancelledError()
                        raise
                    results[key] = result
            finally:
                if cancel_token is not None:
                    cancel_token.clear_pool()

    return results
