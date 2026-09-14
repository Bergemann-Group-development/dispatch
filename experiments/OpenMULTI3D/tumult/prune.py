"""Per-node line-list pruning through ``hop_prune.x`` child processes.

A synthesis on one atmosphere node reads a few hundred transitions of a master list that
holds hundreds of thousands; which ones is decided by HOP's negligibility test on that
node's thermodynamic states. ``hop_prune.x`` (``HOP/standalone/hop_prune.f90``) makes that
selection once per node from the same namelist groups the runs use and writes it as an input
line list, so every run on the node reads the pruned file as its ``linelist_file``.

This module assembles the pruning namelist from a run's settings (:func:`prune_namelist`),
runs one child per uncached node on the shared :class:`tumult.synthesis.runtime`
execution service (:func:`prune_batch`) and enforces the process contract: a child owns
one HOP instance and one CPU slot, cancellation kills it, and a cancelled or failed child
never publishes its output (the executable writes ``<output>.<pid>.partial`` and renames on
success; the partial is removed here on failure). Each output carries a
``<output>.provenance.nml`` sidecar, readable with :func:`provenance`.

Running HOP in a child rather than in this process removes the process-global HOP
conflict for pruning: :data:`tumult.synthesis.runtime.HOP_LOCK` stays necessary only for
callers that still use HOP in-process, such as the GUI's line-strength bands.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from .namelist import NameList

CHECKOUT = Path(__file__).resolve().parents[1]
EXECUTABLE = Path('HOP') / 'standalone' / 'build' / 'hop_prune.x'
SAMPLE_COUNT = 8      # depths sampled per 1D atmosphere; see hop_prune.f90 &prune_params
PROVENANCE_SUFFIX = '.provenance.nml'


def executable(cwd=None):
    """The ``hop_prune.x`` binary: ``$HOP_PRUNE``, else the standalone build directory
    under ``cwd`` or under this checkout. Raises ``FileNotFoundError`` with the build
    command when none exists; pruning has no in-process fallback."""
    candidates = []
    if os.environ.get('HOP_PRUNE'):
        candidates.append(Path(os.environ['HOP_PRUNE']))
    if cwd is not None:
        candidates.append(Path(cwd) / EXECUTABLE)
    candidates.append(CHECKOUT / EXECUTABLE)
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        f'hop_prune.x not found ({", ".join(str(c) for c in candidates)}); build it with '
        f'"make -C HOP/standalone hop_prune.x" or point $HOP_PRUNE at it')


def prune_namelist(output_file, atmos_file, linelist_files, window_file, *, R, dlam,
                   vmic=None, feh=None, abundances='', relative_abundances=True,
                   base=None, sample_count=SAMPLE_COUNT, atmos_format=None, use_abnd=None):
    """The ``hop_prune.x`` input for one atmosphere node.

    ``base`` is the runs' own settings; its ``composition_params``, ``opac_table_params``
    and ``linelist_params`` are carried over so the selection sees the same EOS, background
    vmic, rejection tolerance and element filters the runs do. The node-specific values
    replace what the runs would set per node: ``atmos_file``, ``feh``, ``vmic``, and the
    ``abundances`` override the caller prunes at (the strongest rung of its ladder; a blank
    string leaves the standard composition). ``linelist_files`` and ``dlam`` replace the
    base list selection, and ``window_file`` with ``R`` define the wavelength grid.
    """
    base = base or {}
    nml = NameList()
    atmos = {'atmos_file': str(atmos_file)}
    for key in ('atmos_format', 'use_rho', 'use_ne', 'use_abnd'):
        if key in (base.get('atmos_params') or {}):
            atmos[key] = base['atmos_params'][key]
    if atmos_format:
        atmos['atmos_format'] = atmos_format
    if feh is not None:
        atmos['feh'] = float(feh)
    if vmic is not None:
        atmos['vmic'] = float(vmic)
    if use_abnd is not None:
        atmos['use_abnd'] = bool(use_abnd)
    nml['atmos_params'] = atmos
    composition = dict(base.get('composition_params') or {})
    composition['abundances'] = str(abundances or '')
    composition['relative_abundances'] = bool(relative_abundances)
    nml['composition_params'] = composition
    if base.get('opac_table_params'):
        nml['opac_table_params'] = dict(base['opac_table_params'])
    linelist = dict(base.get('linelist_params') or {})
    linelist.pop('linelist_folder', None)
    linelist['linelist_file'] = [str(f) for f in linelist_files]
    linelist['dlam'] = float(dlam)
    nml['linelist_params'] = linelist
    spectrum = {'R': float(R), 'window_file': str(window_file)}
    if 'in_air' in (base.get('spectrum_params') or {}):
        spectrum['in_air'] = base['spectrum_params']['in_air']
    nml['spectrum_params'] = spectrum
    nml['prune_params'] = {'output_file': str(output_file), 'sample_count': int(sample_count)}
    return nml


def run_prune(nml_path, output_file, *, cwd=None, cancel_token=None, log_file=None):
    """Run ``hop_prune.x`` on ``nml_path`` to completion.

    The child is its own process group so cancellation kills it whole. A nonzero exit
    raises ``RuntimeError`` carrying the child's last diagnostic line; cancellation raises
    :class:`tumult.synthesis.cancel.CancelledError`. Either way the child's partial output
    is removed, so a later caller never mistakes it for a finished list.
    """
    from tumult.synthesis.cancel import CancelledError, kill_process
    output_file = Path(output_file)
    cwd = Path(cwd or CHECKOUT)
    cmd = [str(executable(cwd)), str(Path(nml_path).resolve())]
    log_file = Path(log_file) if log_file else output_file.with_suffix('.log')
    log_file.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, 'OMP_NUM_THREADS': '1'}      # one HOP instance, one CPU slot
    if cancel_token is not None:
        cancel_token.check()
    with open(log_file, 'w') as log:
        process = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log, stderr=log,
                                   start_new_session=(os.name == 'posix'))
        detach = cancel_token.watch_process(process) if cancel_token else (lambda: None)
        try:
            process.wait()
            if cancel_token is not None:
                cancel_token.check()
            if process.returncode:
                raise RuntimeError(f'hop_prune.x exited with status {process.returncode}: '
                                   f'{_last_line(log_file)} ({log_file})')
        except BaseException as exc:
            if not isinstance(exc, CancelledError):
                kill_process(process)
            process.wait()
            partial = output_file.with_name(f'{output_file.name}.{process.pid}.partial')
            partial.unlink(missing_ok=True)
            raise
        finally:
            detach()
    return output_file


def _last_line(log_file):
    try:
        lines = Path(log_file).read_text(errors='replace').strip().splitlines()
    except OSError:
        return ''
    return lines[-1] if lines else ''


def prune_batch(jobs, *, cwd=None, cancel_token=None, runtime=None):
    """Prune every ``(output_file, namelist)`` job whose output does not exist yet.

    Each job becomes one ``hop_prune.x`` child on the execution service, keyed by its
    output path so two workflows asking for the same node share one process, and the
    service's CPU budget bounds how many run at once. Blocks until all are done; the first
    failure cancels the rest. Returns the output paths, cached ones included.
    """
    from tumult.synthesis.cancel import CancelToken
    from tumult.synthesis.runtime import get_runtime

    outputs = []
    pending = []
    for output_file, nml in jobs:
        output_file = Path(output_file)
        outputs.append(output_file)
        if output_file.is_file():
            continue
        output_file.parent.mkdir(parents=True, exist_ok=True)
        nml_path = output_file.with_suffix('.nml')
        # a trailing newline: gfortran does not see a final group whose '/' ends the file
        nml_path.write_text(NameList(nml).to_string() + '\n')
        pending.append((output_file, nml_path, json.dumps(nml, sort_keys=True, default=str)))
    if not pending:
        return outputs

    runtime = runtime or get_runtime()
    batch_token = CancelToken()
    detach = cancel_token.add_callback(batch_token.cancel) if cancel_token else (lambda: None)
    owner = cancel_token or batch_token
    futures = []
    try:
        for output_file, nml_path, signature in pending:
            def run(token, slots, nml_path=nml_path, output_file=output_file):
                return run_prune(nml_path, output_file, cwd=cwd, cancel_token=token)
            futures.append(runtime.submit(str(output_file.resolve()), signature, owner, run,
                                          threads=1, token=batch_token))
        for future in futures:
            future.result()
    finally:
        batch_token.cancel()
        detach()
    return outputs


def provenance(output_file):
    """The provenance record ``hop_prune.x`` wrote next to ``output_file``: the input
    groups as it read them, and ``prune_result`` with the sampling, the effective vmic and
    [Fe/H], the vacuum windows and the number of transitions written."""
    return NameList.from_file(str(output_file) + PROVENANCE_SUFFIX)
