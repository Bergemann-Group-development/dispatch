import os
import copy
import subprocess as subproc
from .auxfuncs import check_file, update_dict
from .tumult import tumult
import re
import numpy as np
import pathlib
from tumult.NameList import NameList

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=None,
           log_file=None, slurm_script=None, verbose=2, dry=False, subprocess=False, runner=None):
    
    input_name = pathlib.Path(input_name)
    input_name.parent.mkdir(parents=True, exist_ok=True)

    cwd = os.getcwd()
    check_file(cwd + '/dispatch.x')
    input = NameList()

    if default is not None:
        if not isinstance(default, NameList) and isinstance(default, (str, pathlib.Path)):
            default = pathlib.Path(default)
            params_nml = default / 'params.nml' if default.is_dir() else default
            default = NameList.from_file(params_nml)
            
        update_dict(input, default)
    
    if nml_settings is not None:
        update_dict(input, nml_settings)

    datadir = input.get('io_params', {'datadir':'data'}).get('datadir', 'data').strip()
    output_folder = datadir / input_name

    if log_file is None and verbose < 2:
        log_file = '/dev/null 2>&1'

    ok = overwrite
    if not ok:
        try:
            run = tumult.read(output_folder, read_atmos=False, read_atom=False)
            existed = True
        except:
            ok = True

    if ok:
        input_file = input_name.with_suffix('.nml')
        input.to_file(input_file)

        if not dry:
            if slurm_script is None:
                if log_file is None:
                    command = f'./dispatch.x {input_file}'
                else:
                    command = f'./dispatch.x {input_file} > {log_file}'
            else:
                command = f'sbatch {slurm_script} {input_file}'

            if verbose == 1:
                print(input_file, ' launched')

            os.environ['OMP_NUM_THREADS'] = str(threads) if threads else ''
            os.environ['OMP_STACKSIZE'] = '256M'
            env = os.environ.copy()

            if runner:
                process = runner.run_command(command, shell=True, env=env, text=True)
                return output_folder, process

            elif subprocess:
                process = subproc.Popen(command, shell=True, env=env)
                return output_folder, process

            else:
                iostat = os.system(command)

            if iostat!=0 and log_file is not None: 
                with open(log_file) as f:
                    for line in f:
                        pass
                    last_line = line
                    print(last_line)

        elif verbose > 0:
            print(input_file, ' written')
            
    else:
        if verbose > 0:
            print(output_folder, ' already exists')

    return output_folder, ok or existed


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