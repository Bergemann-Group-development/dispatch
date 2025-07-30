import os
import f90nml
import copy
import subprocess
from .auxfuncs import check_file
from .m3dis import m3dis

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=None,
           log_file=None, slurm_script=None, verbose=2, dry=False, subprocess=False):
    cwd = os.getcwd()
    check_file(cwd + '/dispatch.x')
    if default is not None:
        if type(default) is f90nml.namelist.Namelist:
            input = copy.deepcopy(default)
        elif default.endswith('.nml'):
            input = f90nml.read(default)
        else:
            input = f90nml.read(default + '/params.nml')
    else:
        input = f90nml.Namelist()
    
    if nml_settings is not None:
        input.patch(nml_settings)

    datadir = 'data/'
    if 'io_params' in input:
        if 'datadir' in input['io_params']:
            datadir = input['io_params']['datadir']

    input_file = input_name + '.nml'
    input_folder = os.path.dirname(input_file)
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)
    output_folder = os.path.join(datadir.strip(), input_name)

    ok = overwrite
    if not ok:
        try:
            run = m3dis.read(output_folder, read_atmos=False, read_atom=False)
            existed = True
        except:
            ok = True

    if ok:
        input.write(input_file, force=True)

        if not dry:
            if slurm_script is None:
                if log_file is None:
                    command = './dispatch.x {:}'.format(input_file)
                else:
                    command = './dispatch.x {:} > {:}'.format(input_file, log_file)
            else:
                command = 'sbatch {:} {:}'.format(slurm_script, input_file)

            if verbose == 1:
                print(input_file, ' launched')

            if verbose < 2 and log_file is None:
                command = command + ' > /dev/null 2>&1'

            sthreads = str(threads) if threads else ''
            os.environ['OMP_NUM_THREADS'] = sthreads
            os.environ['OMP_STACKSIZE'] = '256M'

            if subprocess:
                env = os.environ.copy()
                process = subprocess.Popen(command, shell=True, env=env)
                iostat = process.returncode
            else:
                iostat = os.system(command)

            if iostat!=0 and log_file is not None: 
                with open(log_file) as f:
                    for line in f:
                        pass
                    last_line = line
                    print(last_line)
        elif verbose>0:
            print(input_file, ' written')
    else:
        if verbose>0:
            print(output_folder, ' already exists')

    return output_folder, ok or existed


def monitor(folder):
    n = 0
    nmax = 1

    file = folder + os.sep + 'rank_00000.que'
    if os.path.isfile(file):
        f = open(file, mode='r')
        lines = f.readlines()
        lines = [l for l in lines if 'nlte_patch' in l and '-' in l]
        nmax = len(lines)

    file = folder + os.sep  + 'rank_00000.disp'
    if os.path.isfile(file):
        f = open(file, mode='r')
        lines = f.readlines()
        lines = [l for l in lines if 'END   nlte_patch' in l]
        n = len(lines)

    return n, nmax