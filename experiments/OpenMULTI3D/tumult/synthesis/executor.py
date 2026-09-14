"""Synthesis batches on the shared, budgeted dispatch execution service."""
from __future__ import annotations

import copy
from collections import deque
from concurrent.futures import wait, FIRST_COMPLETED
import json
from pathlib import Path

import tumult as tm
from tumult.auxfuncs import update_dict
from tumult.launch import run_complete, prepare_input

from .cancel import CancelToken
from .runtime import get_runtime, note_run


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
    """Submit every job of a batch with a balanced share of the free CPU budget.

    Follows the ready-queue allocation of ``plan/multi-processing.md``: with ``F``
    slots free (here, this batch's own ceiling ``min(max_threads, runtime.budget)``)
    and ``n`` jobs ready to launch, ``width = min(n, F)`` of them start at once, ``F``
    slots split as evenly as possible (``base = F // width``, the first ``F % width``
    getting one extra). A job's allocation is fixed once it starts; as jobs finish and
    free their slots, the same split is recomputed for whatever is still waiting. Two
    jobs sharing one output (the same node requested by two callers) are always
    submitted as separate subscribers of that shared task — see
    :class:`runtime.ExecutionService` — so only the first occurrence of an output
    consumes a slot from this batch's own accounting; duplicates ride along free.
    Jobs whose output is already a complete, reusable run (and every caller asking for
    it) are submitted at once outside the balanced split, since reading a finished run
    costs no native CPU slot.
    """
    if not jobs:
        return {}
    runtime = runtime or get_runtime()
    batch_token = CancelToken()
    detach = cancel_token.add_callback(batch_token.cancel) if cancel_token else lambda: None
    max_threads = max(1, int(max_threads))
    owner = cancel_token or batch_token
    pending = {}         # future -> job key
    results = {}
    by_output = {}        # output -> [(key, name, nml, kwargs, signature), ...]
    cached_outputs = []
    waiting_outputs = deque()
    future_output = {}    # future -> output, to release that output's slots exactly once

    def submit_entry(name, nml, kwargs, signature, output, slots):
        def run(token, slots=slots):
            token.check()
            options = dict(kwargs, threads=slots)
            # External runners retain their own execution/cancellation contract.
            # Hard process cancellation is available for local synchronous jobs.
            if not any(options.get(k) for k in ('runner', 'slurm_script', 'run_async')):
                options['cancel_token'] = token
            result = _run_single(None, name, nml, options)[1]
            token.check()
            return result
        return runtime.submit(output, signature, owner, run, threads=slots, token=batch_token)

    def submit_output(output, slots):
        for key, name, nml, kwargs, signature in by_output[output]:
            future = submit_entry(name, nml, kwargs, signature, output, slots)
            pending[future] = key
            future_output[future] = output

    free = max(1, min(max_threads, runtime.budget))
    launched_slots = {}    # output -> slots granted from `free` (cache hits are not here)

    def allocate_ready():
        nonlocal free
        n = len(waiting_outputs)
        if n == 0 or free <= 0:
            return
        width = min(n, free)
        base, remainder = divmod(free, width)
        for i in range(width):
            output = waiting_outputs.popleft()
            slots = base + (1 if i < remainder else 0)
            launched_slots[output] = slots
            free -= slots
            submit_output(output, slots)

    try:
        # Resolve every job's output/signature, own its input snapshot (workflows may
        # keep changing their live settings after this), and group by output so a node
        # requested twice in one batch shares one native run without inflating the
        # count of "ready" jobs the width/base/remainder split below is computed
        # against. Inside the try so a cancellation here still detaches cleanly.
        for key, name, nml in jobs:
            batch_token.check()
            nml = copy.deepcopy(nml)
            kwargs = copy.deepcopy(launch_kwargs)
            prepared = prepare_input(nml, kwargs.get('default'))
            output = str((Path(kwargs.get('cwd') or '.') /
                          _output_folder(name, prepared)).resolve())
            # CPU allocation is scheduling metadata, not part of a spectrum's identity.
            signature = json.dumps([prepared, {k: v for k, v in kwargs.items()
                                               if k not in ('threads', 'verbose', 'log_file')}],
                                   sort_keys=True, default=str)
            entry = (key, name, nml, kwargs, signature)
            if output not in by_output:
                by_output[output] = [entry]
                note_run(output)
                if not kwargs.get('overwrite') and run_complete(output):
                    cached_outputs.append(output)
                else:
                    waiting_outputs.append(output)
            else:
                by_output[output].append(entry)

        for output in cached_outputs:
            submit_output(output, 1)      # a cache read costs no native CPU slot
        allocate_ready()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            batch_token.check()
            for future in done:
                key = pending.pop(future)
                results[key] = future.result()
                output = future_output.pop(future)
                slots = launched_slots.pop(output, None)  # None for a cache-hit output
                if slots is not None:
                    free += slots
            allocate_ready()
        return results
    finally:
        batch_token.cancel()  # also detach outstanding subscribers on an error
        detach()
