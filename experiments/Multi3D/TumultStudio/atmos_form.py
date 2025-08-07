from bokeh.models import *
from . import ServerManager
from .ServerManager import update_data_source
from .auxfuncs import *
from .StudioAtmos import StudioAtmos
import numpy as np
import ipywidgets as widgets
from .settings import saved_settings
import ipyformkit as ifk
from .synth import run_synth
import pathlib

fields = {'atmosphere format':('Marcs', 'Must', 'Stagger', 'Stagger2', 'Text', 'Multi1D', 'Siu_bin'),
          'atmosphere file': ('no atmosphere picked',),
          ('nx', 'ny'): (10, 10),
          ('FeH', 'vmic'): (0.0, -1.0),
          'import':None}

tooltips = {
    'FeH':'Metallicity of the atmosphere, zero being solar metallicity',
    'vmic':'Mircoturbulence velocity in km/s. Default is -1, which will adapt vmic if given by the atmosphere file and if not it is zero.',
    }

mhd_formats = ('Must', 'Stagger', 'Stagger2')
hide = {'nx': lambda d: d['atmosphere format'] not in mhd_formats,
        'ny': lambda d: d['atmosphere format'] not in mhd_formats,}

mandatory = ['atmosphere file']
atmos_import_form = ifk.Form(fields, mandatory=mandatory, tooltips=tooltips, hide=hide)
import_button = atmos_import_form.widgets_dict['import'].wid
atmos_wid = atmos_import_form.widgets_dict['atmosphere file'].wid
atmos_fmt = atmos_import_form.widgets_dict['atmosphere format'].wid

@log_callback_errors
def list_atmospheres(change, wid):
    fmt = change['new']
    folder = pathlib.Path(saved_settings['atmospheres folder'])
    
    # Get all files in the folder
    files = list(folder.iterdir())
    
    if fmt == 'Marcs':
        files = [f for f in files if f.suffix == '.mod']
    elif fmt == 'Multi1D':
        files = [f for f in files if f.name.startswith('atmos.')]
    elif fmt == 'Siu_bin':
        files = [f for f in files if f.suffix == '.dat']
    elif fmt == 'Must':
        files = get_common_suffix_files(files, ['.atmos', '.mesh'])
    
    # Convert to relative paths (names only)
    files = [f.name for f in files]
    
    wid.options = files
    wid.value = files[0] if files else None

atmos_fmt.observe(lambda change: list_atmospheres(change, atmos_wid), names='value')
list_atmospheres({'new':atmos_fmt.value}, atmos_wid)

def init_import_atmos_button(studio):
    import_atmos_button = widgets.Button(description='Import Atmos')
        
    @log_callback_errors
    def on_import_atmos_click(btn):
        studio.browser.selected = None
        studio.right_sidebar.children = (atmos_import_form.vbox,)
        studio.right_sidebar.titles = ['Atmos Settings']

    @log_callback_errors
    def on_import_click(btn):
        import_atmosphere(studio)
    
    import_atmos_button.on_click(on_import_atmos_click)
    import_button.on_click(on_import_click)
    studio.controls.children += (import_atmos_button,)

def get_atmos_settings():
    values = atmos_import_form.check_and_return_values()
    atmos_folder = pathlib.Path(saved_settings['atmospheres folder'])
    settings = {'atmos_params':{'atmos_file':atmos_folder / values['atmosphere file'],
                                'atmos_format':values['atmosphere format'],
                                'feh':values['FeH'],
                                'vmic':values['vmic'],
                               },
                'm3d_params':{'maxiter':0,},
                'linelist_params':{'line_lists':''},
                }
    
    if values['atmosphere format'] in mhd_formats:
        settings['atmos_params']['nx'] = values['nx']
        settings['atmos_params']['ny'] = values['ny']

    return settings

def import_atmosphere(studio):
    atmospheres = studio._get_objects_by_type(StudioAtmos)
    i = find_smallest_missing([a.ic for a in atmospheres])
    name = f'Atmos_{i:03d}'
    settings = get_atmos_settings()
    folder = studio.synth_folder
    studio.progress.description = f'Importing atmosphere: {name}...'
    run = run_synth(folder, name, settings, progress=studio.progress)
    new_atmos = StudioAtmos.from_run(run, studio=studio)
    new_atmos.ic = i
    studio.browser.add_object(new_atmos)
    studio.progress.description = f'Done!'