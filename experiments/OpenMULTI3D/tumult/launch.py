import os
import copy
import subprocess
from .auxfuncs import check_file, update_dict
import re
import numpy as np
import pathlib
from tumult.namelist import NameList

def _load_default(item):
    """Load one default into a NameList-compatible mapping.

    ``item`` may be a NameList, a path to a .nml file, or a directory containing
    ``params.nml``.
    """
    if isinstance(item, NameList):
        return item
    p = pathlib.Path(item)
    params = p / 'params.nml' if p.is_dir() else p
    return NameList.from_file(params)

def prepare_input(nml_settings, default):
    """Handles logic for creating the NameList configuration.

    ``default`` may be a single default (NameList, .nml path, or directory) or a
    list/tuple of them, applied in order so each later entry overrides earlier ones.
    """
    input_data = NameList()

    # Handle default settings
    if default is not None:
        items = default if isinstance(default, (list, tuple)) else [default]
        for item in items:
            update_dict(input_data, _load_default(item))

    # Apply specific overrides
    if nml_settings is not None:
        update_dict(input_data, nml_settings)

    return input_data

def run_complete(folder):
    """Whether ``folder`` holds a finished, reusable dispatch run.

    dispatch.x creates its run directory (and the ``.log`` / ``rank_*`` files) up front and
    only writes ``abund.txt`` once composition setup has succeeded, so a folder that exists
    but lacks ``abund.txt`` is a crashed or half-written run. Treating such a folder as a cache
    hit both crashes the reader (``tumult.read`` needs ``abund.txt``) and poisons the cache —
    every later synthesis that selects that node reuses the broken folder until the cache is
    cleared. So the cache must key on this predicate, not on directory existence alone.
    """
    return (pathlib.Path(folder) / 'abund.txt').exists()


def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=None, cwd=None,
           log_file=None, slurm_script=None, verbose=2, dry=False, run_async=False, runner=None,
           cancel_token=None):

    if cancel_token is not None:
        cancel_token.check()
        if runner or run_async or slurm_script:
            raise ValueError("Managed cancellation requires a synchronous local dispatch run")
    
    # Path setup
    input_path = pathlib.Path(input_name)
    cwd = pathlib.Path(cwd or pathlib.Path.cwd()).expanduser().resolve()
    (cwd / input_path).parent.mkdir(parents=True, exist_ok=True)
    executable = cwd / 'dispatch.x'

    # Build configuration and determine output folder
    nml = prepare_input(nml_settings, default)
    datadir = nml.get('io_params', {}).get('datadir', 'data').strip()
    output_folder = cwd / datadir / input_path

    # Check for existing data. Only a complete run counts as a cache hit; a partial folder
    # from a crashed run is re-run (dispatch overwrites its outputs) rather than served back.
    if output_folder.exists() and run_complete(output_folder) and not overwrite:
        if verbose > 0:
            print(f"{output_folder} already exists")
        return output_folder, False

    # Write the NML file
    input_file = input_path.with_suffix('.nml')
    nml.to_file(cwd / input_file)
    
    if dry:
        if verbose > 0: print(f"{input_file} written (dry run)")
        return output_folder, True

    # Prepare execution environment
    env = {**os.environ, 'OMP_STACKSIZE': '256M'}
    if threads:
        env['OMP_NUM_THREADS'] = str(threads)

    # Determine command and output handling
    if slurm_script:
        cmd = ['sbatch', slurm_script, str(input_file)]
    else:
        cmd = [str(executable), str(input_file)]

    # Open log file if requested
    if verbose < 2 and log_file is None:
      log_file = output_folder.with_suffix('.log')

    if log_file is not None:
        # dispatch.x creates output_folder itself only once it runs (after this point),
        # so the log file's parent directory may not exist yet — create it first.
        log_file = pathlib.Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

    log_handle = open(log_file, 'w') if log_file else None

    if verbose == 1:
        print(f"{input_file} launching")

    # tumult.synthesis.jobrecord reads run_complete from this module, so it is imported
    # here rather than at the top: a run's record of itself belongs to the synthesis
    # layer, and the launcher only needs it while it is starting a run.
    from tumult.synthesis import jobrecord
    from tumult.synthesis.cancel import CancelledError, kill_process

    # Execution logic
    try:
        # An externally scheduled or asynchronous run keeps its own execution contract:
        # this call returns before the run does, so the record stays 'submitted' and the
        # holder of the returned handle owns the ending. A local synchronous run is
        # recorded as running and finished here.
        if runner:
            jobrecord.write(output_folder, nml, state='submitted')
            return output_folder, runner.run_command(cmd, env=env)

        if run_async:
            process = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log_handle, stderr=log_handle)
            jobrecord.write(output_folder, nml, pid=process.pid, state='submitted')
            return output_folder, process

        if cancel_token is None:
            process = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log_handle, stderr=log_handle)
            # sbatch exits as soon as the job is queued, so its pid is not the run's and
            # its exit status does not say whether the run succeeded.
            jobrecord.write(output_folder, nml, pid=None if slurm_script else process.pid,
                            state='submitted' if slurm_script else 'running')
            process.wait()
            result = process
            if not slurm_script:
                jobrecord.finish(output_folder,
                                 'finished' if result.returncode == 0 else 'failed')
        else:
            cancel_token.check()
            process = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log_handle, stderr=log_handle,
                                       start_new_session=(os.name == 'posix'))
            jobrecord.write(output_folder, nml, pid=process.pid)
            detach = cancel_token.watch_process(process)
            try:
                process.wait()
                # The termination cause is read here, before any cleanup touches the token:
                # a cancelled token means the user (or session) ended the run, whatever the
                # killed child's exit status says; otherwise a nonzero status is a failure.
                cancel_token.check()
                if process.returncode:
                    raise RuntimeError(f"dispatch exited with status {process.returncode}: {log_file}")
            except BaseException as exc:
                cancelled = isinstance(exc, CancelledError)
                if not cancelled:
                    kill_process(process)   # an interrupted wait must not leave the child running
                process.wait()       # reap the child before releasing its CPU allocation
                # abund.txt is written during initialisation, before completion.
                # Failed/cancelled jobs must not become cache hits on the next run.
                (output_folder / 'abund.txt').unlink(missing_ok=True)
                jobrecord.finish(output_folder, 'cancelled' if cancelled else 'failed',
                                 returncode=process.returncode,
                                 log=str(log_file) if log_file else None,
                                 error=None if cancelled else f'{type(exc).__name__}: {exc}')
                raise
            finally:
                detach()
            result = process
            jobrecord.finish(output_folder, 'finished')

        if result.returncode != 0 and log_file:
            print(f"Error: Command failed. Check {log_file} for details.")
            
    finally:
        if log_file is not None:
            log_handle.close()

    return output_folder, True


# The per-atom convergence line dispatch.x writes once per iteration:
#   local = -3.86E-05, global = -1.22E-03, lvl =   1, pos=[...] Mg
# A bare `global =` also matches the `flux global = ...` line that follows it, which would
# count every iteration twice; anchoring on `local =` picks the population series, which is
# the one the `Converged` banner is gated on.
GLOBAL_ERROR = re.compile(r'^\s*local\s*=.*?\bglobal\s*=\s*([^,]+)', re.MULTILINE)


def monitor(folder, log_file=None):
    """Progress of a running dispatch job, as ``(done, total)``, from its own output files.

    Two sources, added together. The log's per-iteration convergence errors give the NLTE
    solve's progress against ``maxiter``, replaced by an extrapolated total once the last
    few iterations are improving monotonically enough to fit. The dispatcher's queue and
    dispatch logs give the patch tasks, which for an NLTE run means the atmosphere setup
    only -- the iteration tasks are not written in a form these files distinguish, so they
    contribute a fixed count rather than a rate.

    Everything is optional: a run that has not written its namelist yet, or one launched
    without a log file, reports what it can rather than failing. ``total`` is an estimate
    and ``done`` can exceed it when a solve outruns its own extrapolation, so a caller
    rendering a bar must clamp the ratio.
    """
    folder = pathlib.Path(folder)
    nml_file = folder / 'params.nml'
    n = 0
    nmax = 0

    if nml_file.is_file() and log_file and pathlib.Path(log_file).is_file():
        nml = NameList.from_file(nml_file)
        maxiter = nml['rt_iteration_params'].get('maxiter', 0)
        convlim = nml['nlte_params'].get('convlim', None)
        if convlim is not None and maxiter > 0:
            convlim = np.log10(convlim)
            content = pathlib.Path(log_file).read_text(errors='replace')
            # the first entry is the LTE iteration's global = 0, whose log is -inf
            matches = GLOBAL_ERROR.findall(content)[1:]
            maxerr = convlim - np.log10([abs(float(m)) for m in matches]) if matches \
                else np.empty(0)
            latest = maxerr[-5:]
            nl = len(latest)
            if len(maxerr) >= 3 and all(np.diff(latest) > 0):
                p = np.polyfit(latest, np.arange(nl), deg=2)
                ndone = len(maxerr)
                nleft = np.ceil(np.polyval(p, 0)) - nl + 1
                n += ndone
                nmax += ndone + nleft
            else:
                n += len(maxerr)
                nmax += maxiter

    # The `-` selects output patches, which carry negative frequency indices (see
    # docs/GLOSSARY.md). An NLTE solve's own nlte_patch lines carry none, so on that run
    # shape these two counts are the eight mhd_patch setup tasks and nothing else.
    file = folder / 'rank_00000.que'
    if file.is_file():
        lines = file.read_text(errors='replace').splitlines()
        nmax += len([l for l in lines if 'mhd_patch'  in l])
        nmax += len([l for l in lines if 'nlte_patch' in l and '-' in l])

    file = folder / 'rank_00000.disp'
    if file.is_file():
        lines = file.read_text(errors='replace').splitlines()
        n += len([l for l in lines if 'END    mhd_patch' in l])
        n += len([l for l in lines if 'END   nlte_patch' in l and '-' in l])

    return n, nmax
