"""Synthesis batches on the shared, budgeted dispatch execution service."""
from __future__ import annotations

import copy
from concurrent.futures import wait, FIRST_COMPLETED
import json
from pathlib import Path

import tumult as tm
from tumult.auxfuncs import update_dict
from tumult.launch import run_complete, prepare_input

from .cancel import CancelToken
from .runtime import get_runtime


def _run_single(key, nml_name, nml, launch_kwargs):
    """Launch/read one job in a supervisor thread."""
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


def _run_batch(jobs, launch_kwargs, max_threads, cancel_token=None, *, runtime=None):
    """Submit a bounded batch; preserve job keys and deduplicate shared outputs.

    The service owns the global CPU budget. ``max_threads`` additionally limits
    this batch's outstanding work. Cached outputs are read by the same supervisor
    path, so a file written by an in-flight job is never mistaken for a cache hit.
    """
    if not jobs:
        return {}
    runtime = runtime or get_runtime()
    batch_token = CancelToken()
    detach = cancel_token.add_callback(batch_token.cancel) if cancel_token else lambda: None
    max_threads = max(1, int(max_threads))
    threads = min(max_threads, max(1, int(launch_kwargs.get("threads") or 1)))
    outstanding = max(1, min(max_threads, runtime.budget) // threads)
    owner = cancel_token or batch_token
    pending = {}
    remaining = iter(jobs)
    results = {}

    def submit_next():
        batch_token.check()
        try:
            key, name, nml = next(remaining)
        except StopIteration:
            return False
        # Own the input snapshot: workflows may continue changing their settings.
        nml = copy.deepcopy(nml)
        kwargs = copy.deepcopy(launch_kwargs)
        prepared = prepare_input(nml, kwargs.get('default'))
        output = str(_output_folder(name, prepared).resolve())
        # CPU allocation is scheduling metadata, not part of a spectrum's identity.
        signature = json.dumps([prepared, {k: v for k, v in kwargs.items()
                                           if k not in ('threads', 'verbose', 'log_file')}],
                               sort_keys=True, default=str)

        def run(token, slots):
            token.check()
            options = dict(kwargs, threads=slots)
            # External runners retain their own execution/cancellation contract.
            # Hard process cancellation is available for local synchronous jobs.
            if not any(options.get(k) for k in ('runner', 'slurm_script', 'run_async')):
                options['cancel_token'] = token
            result = _run_single(None, name, nml, options)[1]
            token.check()
            return result

        future = runtime.submit(output, signature, owner, run,
                                threads=threads, token=batch_token)
        pending[future] = key
        return True

    try:
        for _ in range(outstanding):
            if not submit_next():
                break
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            batch_token.check()
            for future in done:
                key = pending.pop(future)
                results[key] = future.result()
                submit_next()
        return results
    finally:
        batch_token.cancel()  # also detach outstanding subscribers on an error
        detach()
