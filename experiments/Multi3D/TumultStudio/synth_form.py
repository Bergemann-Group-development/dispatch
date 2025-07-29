from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .LineMask import LineMask
from .llist_annotations import get_llist_table
from m3dis.auxfuncs import update_dict
import ipyformkit as ifk
import m3dis
import numpy as np
import tempfile
import glob

fields = {('line_mask', 'atmosphere format'):(('no line masks yet',),
                                             ('Marcs', 'Must', 'Stagger', 'Stagger2', 'Text', 'Multi1D', 'Siu_bin')),
          'atmosphere file': ('no atmosphere picked',),
          ('FeH', 'vmic'): (0.0, -1.0),
          ('abundances', 'relative'):('O=8.70, Fe=7.50', False),
          ('exclude elements','only elements'):('O, Fe', 'O, Fe'),
          'launch':None}

tooltips = {
    'FeH':'Metallicity of the atmosphere, zero being solar metallicity',
    'vmic':'Mircoturbulence velocity in km/s. Default is -1, which will adapt vmic if given by the atmosphere file and if not it is zero.',
    'relative':'If true please specify abundances as [X/Fe], otherwise please specify absolute abundances X/H',
    'exclude elements':'Elements to exclude from the linelist.',
    'only elements':'All other elements will be removed from the linelist. If empty, all elements will be included.',
    }

mandatory = ['atmosphere file']

check = {'line_mask': lambda d: type(d['line_mask'])==LineMask}

synth_form = ifk.Form(fields, mandatory=mandatory, check=check, tooltips=tooltips)

disable = {
    'launch': lambda d: synth_form.check_and_return_values(verbose=False) is None,
    'exclude elements': lambda d: d['only elements'] != '',
    }

synth_form._disable_conditions = synth_form.add_observer(disable, synth_form.update_disable)
synth_form.update_disable()

atmos_wid = synth_form.widgets_dict['atmosphere file'].wid
atmos_fmt = synth_form.widgets_dict['atmosphere format'].wid
launch_btn = synth_form.widgets_dict['launch'].wid

def list_atmospheres(change):
    fmt = change['new']
    folder = saved_settings['atmospheres folder'] + os.sep
    files = glob.glob(folder + '*')
    if fmt == 'Marcs':
        files = [f for f in files if f.endswith('.mod')]
    elif fmt == 'Multi1D':
        files = [f for f in files if f.split(os.sep)[-1].startswith('atmos.')]
    elif fmt == 'Siu_bin':
        files = [f for f in files if f.endswith('.dat')]
    elif fmt == 'Must':
        files = get_common_suffix_files(files, ['.atmos', '.mesh'])

    files = [f[len(folder)-1:] for f in files]
    
    atmos_wid.options = files
    atmos_wid.value = files[0]

atmos_fmt.observe(list_atmospheres, names='value')
list_atmospheres({'new':atmos_fmt.value})

def get_synth_settings(values):
    atmos_file = saved_settings['atmospheres folder'] + os.sep + values['atmosphere file']
    settings = {'atmos_params':{'atmos_file':atmos_file,
                                'atmos_format':values['atmosphere format'],
                                'feh':values['FeH'],
                                'vmic':values['vmic'],
                               },
                'composition_params':{'abundances':values['abundances'],
                                      'relative_abundances':values['relative']},
                'linelist_params':{'exclude_elements':values['exclude elements'],
                                   'only_elements':values['only elements']}}
    return settings

def launch_synth(name, progress=None):
    values = synth_form.check_and_return_values()
    line_mask = values['line_mask']
    settings = get_synth_settings(values)
    run = run_synth(name, settings, line_mask, progress=progress)
    return run

def run_synth(name, settings, line_mask, default=None, llist=None, progress=None):
    folder = saved_settings['synthesis folder'] + os.sep
    os.path.exists(folder) or os.makedirs(folder)
    input_name = folder + name
    lam = line_mask.get_wavelength_grid(daa=0.01)

    if default:
        update_dict(settings, {
            'atmos_params': {'FeH': 0},
            'composition_params': {'abund_file': default + 'abund.txt'}})
    else:
        default = module_dir + 'default.nml'

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as llist_file:
        if llist is not None: llist.write_to_ts_format(llist_file.name)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as lam_file:
            np.savetxt(lam_file.name, lam)
            settings['spectrum_params'] = {'lam_file':lam_file.name}
            if llist is not None: settings['linelist_params'] = {'line_lists':llist_file.name}

            @run_in_thread
            def execute_m3dis():
                folder, ok = m3dis.launch(input_name, nml_settings=settings, default=default,
                                        log_file='test.log', threads=2, overwrite=True, verbose=0)
                return folder, ok
            
            thread = execute_m3dis()

            if progress:
                while thread.is_alive():
                    monitor_synth(name, progress)
                    time.sleep(0.25)

            folder, ok = thread.out['result']
        
    run = m3dis.read(folder, read_llist=True)
    run.line_mask = line_mask
    run.llist_html = get_llist_table(run)
    run.label = name
    return run

def monitor_synth(name, progress):
    folder = 'data' + os.sep + saved_settings['synthesis folder'] + os.sep + name
    n, nmax = m3dis.monitor(folder)
    progress.value = n
    progress.max = nmax+1