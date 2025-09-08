from bokeh.models import *
from . import ServerManager
from types import SimpleNamespace
from .StudioGrid import StudioGrid
from .auxfuncs import *
from .StudioAtmos import StudioAtmos
import numpy as np
import ipywidgets as widgets
from .settings import saved_settings
import ipyformkit as ifk
from .synth import *
import pathlib
import parse

fields = {'atmosphere format':('Marcs', 'Text', 'Multi1D', 'Siu_bin'),
            'grid folder': './input_multi3d/atmos/',
            'pattern preset': (
                '{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}',
                'custom'
                ),
            'pattern':'type parse pattern here', 
            'import grid': void}

hide = {'pattern': lambda d: d['pattern preset'] != 'custom'}
atmos_grid_form = ifk.Form(fields, hide=hide)
import_btn = atmos_grid_form.inputs_dict['import grid']
atmos_fmt = atmos_grid_form.inputs_dict['atmosphere format']
grid_folder = atmos_grid_form.inputs_dict['grid folder']
pattern_wid = atmos_grid_form.inputs_dict['pattern']
preset_wid = atmos_grid_form.inputs_dict['pattern preset']

@log_callback_errors
def on_preset_change(change):
    preset = change['new']
    if preset != 'custom':
        pattern_wid.value = preset
    else:
        pattern_wid.value = 'type parse pattern here'

preset_wid.observe(on_preset_change, names='value')
on_preset_change({'new':preset_wid.value})

def init_interp_atmos_form(studio):
    grid_folder.value = saved_settings['atmosphere grid folder']

    @log_callback_errors
    def on_init_click(btn):
        studio.add_object('AtmosGrid', import_atmosphere_grid)

    import_btn.on_click(on_init_click)

def import_atmosphere_grid(studio, name, runner=None):
    values = atmos_grid_form.check_and_return_values()
    fmt = values['atmosphere format']
    atmos_files = filter_atmospheres(values['grid folder'], fmt)
    p = parse.compile(pattern_wid.value)

    data = []
    files = []
    for f in atmos_files:
        result = p.search(f.name)
        if result:
            files.append(f)
            data.append(result.named)
    
    if len(data) > 0:
        grid = StudioGrid(name, data, files, fmt)
        return grid
    else:
        studio.progress.description = 'No atmospheres found matching the pattern.'
