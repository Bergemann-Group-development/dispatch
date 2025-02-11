import os
import f90nml
import copy
from .auxfuncs import check_file
from .m3dis import m3dis

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=32,
           log_file=None, slurm_script=None, verbose=2, dry=False):
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
    output_folder = os.path.join(datadir.strip(), input_name)

    ok = overwrite
    if not ok:
        try:
            run = m3dis.read(output_folder, read_atmos=False, read_atom=False)
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

            if verbose>1:
                print(input_file, ' launched')
            else:
                command = command + ' > /dev/null 2>&1'

            os.environ['OMP_NUM_THREADS'] = str(threads)
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

    return output_folder, ok