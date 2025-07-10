from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .LineMask import LineMask
import ipyformkit as ifk
import m3dis
import numpy as np
import tempfile
import glob

fields = {('line_mask', 'atmosphere format'):(('no line masks yet',),
                                             ('Marcs', 'Must', 'Stagger', 'Stagger2', 'Text')),
          'atmosphere file': ('no atmosphere picked',),
          'abundances':'O=8.70, Fe=7.50',
          ('exclude elements','only elements'):('O, Fe', 'O, Fe'),
          'launch':None}

mandatory = ['atmosphere file']
check = {'line_mask': lambda d: type(d['line_mask'])==LineMask}
synth_form = ifk.Form(fields, title='Synthesis Settings', mandatory=mandatory, check=check)
atmos_wid = synth_form.widgets_dict['atmosphere file'].wid
atmos_fmt = synth_form.widgets_dict['atmosphere format'].wid

def list_atmospheres(change):
    fmt = change['new']
    folder = saved_settings['atmospheres folder'] + os.sep
    files = glob.glob(folder + '*')
    if fmt == 'Marcs':
        files = [f for f in files if f.endswith('.mod')]
    elif fmt == 'Must':
        files = get_common_suffix_files(files, ['.atmos', '.mesh'])

    files = [f[len(folder)-1:] for f in files]
    
    atmos_wid.options = files
    atmos_wid.value = files[0]

atmos_fmt.observe(list_atmospheres, names='value')
list_atmospheres({'new':atmos_fmt.value})

def get_synth_settings(values, lam_file):
    atmos_file = saved_settings['atmospheres folder'] + os.sep + values['atmosphere file']
    settings = {'atmos_params':{'atmos_file':atmos_file,
                                'atmos_format':values['atmosphere format']},
                'spectrum_params':{'lam_file':lam_file.name},
                'composition_params':{'abundances':values['abundances']},
                'linelist_params':{'exclude_elements':values['exclude elements'],
                                   'only_elements':values['only elements']}}
    return settings

@log_callback_errors
def launch_synth():
    default = 'TumultStudio/default.nml'
    values = synth_form.check_and_return_values()
    line_mask = values['line_mask']
    lam = line_mask.source.data['x'][line_mask.islc]
    nthreads = 2# if values['atmosphere format'] in ['Marcs', 'Text'] else 4
    
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as lam_file:
        np.savetxt(lam_file.name, lam)
        settings = get_synth_settings(values, lam_file)
        folder, ok = m3dis.launch('test', nml_settings=settings, default=default,
                                  log_file='test.log', threads=nthreads, overwrite=True, verbose=0)
    new_run = m3dis.read(folder)
    new_run.line_mask = line_mask
    return new_run
