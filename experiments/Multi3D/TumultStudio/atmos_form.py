from bokeh.models import *
from .auxfuncs import *
from .StudioAtmos import StudioAtmos
import numpy as np
import ipywidgets as widgets
from .settings import saved_settings
import ipyformkit as ifk
from .synth import *
import pathlib
from .atmos_grid_form import atmos_grid_form
from types import SimpleNamespace

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
def update_atmos_fmt(change):
    fmt = change['new']
    folder = saved_settings['atmospheres folder']
    files = filter_atmospheres(folder, fmt)
    # Convert to relative paths (names only)
    atmos_wid.options = [f.name for f in files]
    atmos_wid.value = atmos_wid.options[0] if files else None

atmos_fmt.observe(update_atmos_fmt, names='value')
update_atmos_fmt({'new':atmos_fmt.value})

def init_import_atmos_button(studio):
    import_atmos_button = widgets.Button(description='Import Atmos')
        
    @log_callback_errors
    def on_import_atmos_click(btn):
        studio.browser.selected = None
        studio.right_sidebar.children = (atmos_import_form.vbox, atmos_grid_form.vbox)
        studio.right_sidebar.titles = ['Direct Import', 'Import Grid']

    @log_callback_errors
    def on_import_click(btn):
        studio.add_object('Atmos', import_atmosphere)
    
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

def import_atmosphere(studio, name, runner=None):
    settings = get_atmos_settings()
    run = run_synth(studio, name, settings, runner=runner)
    if run:
        run = StudioAtmos.from_run(run)
        return run
    