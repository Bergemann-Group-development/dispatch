import os
import copy
import subprocess
from .auxfuncs import check_file, update_dict
from .tumult import tumult
import re
import numpy as np
import pathlib
from tumult.NameList import NameList

def prepare_input(nml_settings, default):
    """Handles logic for creating the NameList configuration."""
    input_data = NameList()
    
    # Handle default settings
    if default is not None:
        if not isinstance(default, NameList):
            default_path = pathlib.Path(default)
            params = default_path / 'params.nml' if default_path.is_dir() else default_path
            default = NameList.from_file(params)
        update_dict(input_data, default)
    
    # Apply specific overrides
    if nml_settings is not None:
        update_dict(input_data, nml_settings)
        
    return input_data

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=None, cwd=None,
           log_file=None, slurm_script=None, verbose=2, dry=False, run_async=False, runner=None):
    
    # Path setup
    input_path = pathlib.Path(input_name)
    input_path.parent.mkdir(parents=True, exist_ok=True)
    cwd = pathlib.Path(cwd or os.getcwd())
    executable = cwd / 'dispatch.x'

    # Build configuration and determine output folder
    nml = prepare_input(nml_settings, default)
    datadir = nml.get('io_params', {}).get('datadir', 'data').strip()
    output_folder = pathlib.Path(datadir) / input_path

    # Check for existing data
    if output_folder.exists() and not overwrite:
        if verbose > 0:
            print(f"{output_folder} already exists")
        return output_folder, False

    # Write the NML file
    input_file = input_path.with_suffix('.nml')
    nml.to_file(input_file)
    
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

    # Execution logic
    try:
        if runner:
            return output_folder, runner.run_command(cmd, env=env)

        if run_async:
            return output_folder, subprocess.Popen(cmd, env=env, stdout=log_handle, stderr=log_handle)

        # Standard blocking execution
        result = subprocess.run(cmd, env=env, stdout=log_handle, stderr=log_handle)
        
        if result.returncode != 0 and log_file:
            print(f"Error: Command failed. Check {log_file} for details.")
            
    finally:
        if log_file is not None:
            log_handle.close()

    return output_folder, True


def monitor(folder, log_file=None):
    folder = pathlib.Path(folder)
    nml_file = folder / 'params.nml'
    n = 0
    nmax = 0
    ndone = 0
    nleft = 0

    if os.path.isfile(nml_file):
        f = open(log_file, mode='r')
        nml = NameList.from_file(nml_file)
        if 'atom_params' in nml and 'm3d_params' in nml:
            maxiter = nml['m3d_params']['maxiter']
            convlim = np.log10(nml['atom_params']['convlim'])
            if maxiter > 0:
                if os.path.isfile(log_file):
                    content = f.read()
                    # Find all global values using regex
                    global_pattern = r'global\s*=\s*([^,]+)'
                    matches = re.findall(global_pattern, content)
                    maxerr = convlim - np.log10([abs(float(m)) for m in matches[1:]])
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
        f.close()

    file = folder / 'rank_00000.que'
    if file.is_file():
        f = open(file, mode='r')
        lines = f.readlines()
        nmax += len([l for l in lines if 'mhd_patch'  in l])
        nmax += len([l for l in lines if 'nlte_patch' in l and '-' in l])
        f.close()

    file = folder / 'rank_00000.disp'
    if file.is_file():
        f = open(file, mode='r')
        lines = f.readlines()
        n += len([l for l in lines if 'END    mhd_patch' in l])
        n += len([l for l in lines if 'END   nlte_patch' in l and '-' in l])
        f.close()

    return n, nmax