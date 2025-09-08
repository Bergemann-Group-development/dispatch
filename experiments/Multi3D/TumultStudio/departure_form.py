from .ServerManager import update_source_data, clear_data
from .settings import saved_settings
from .synth import run_synth, StudioSynth
from .StudioDepartures import StudioDeparture
from .auxfuncs import *
from . import ServerManager
from .threading import execute_in_thread
import ipyformkit as ifk
import numpy as np
import os
import logging
import pathlib
import tumult

atom_files = [f for f in os.listdir(saved_settings['atom folder']) if not f.startswith('.')]

fields = {'atom file': tuple(atom_files),
          ('abundance', 'relative'): (0.0, True),
          'Compute Departures': void}

departure_form = ifk.Form(fields)
compute_departures_btn = departure_form.inputs_dict['Compute Departures']

def init_departure_form(studio):
    @log_callback_errors
    def on_compute_click(btn):
        name = departure_form['atom file']
        item = studio.browser.selected
        studio.add_object(name, compute_departures, parent=item)
        
    compute_departures_btn.on_click(on_compute_click)
  
def get_departure_settings(atmos):
    values = departure_form.check_and_return_values()
    atom_folder = pathlib.Path(saved_settings['atom folder'])
    if values:
        name = atmos.label + '_' + values['atom file'].replace('.', '_')
        settings = {'atmos_params':{'atmos_format':'save',
                                    'atmos_file':atmos.folder,
                                },
                    'atom_params':{'atom_file':atom_folder / values['atom file'],
                                    'abundance':values['abundance'],
                                    'relative_abundance':values['relative'],
                                    },
                    }
        return settings, name
    else:
        return None, None
    
def compute_departures(studio, name, runner=None):
    item = studio.browser.selected
    atmos = item.obj
    default = module_dir / 'default_nlte.nml'
    settings, name = get_departure_settings(atmos)
    if settings is None: return

    run = run_synth(studio, name, settings, default=default, runner=runner)
    if run:
        link_files(run.sfolder, atmos.sfolder, 'pop_*')
        link_files(run.folder/'atoms', atmos.folder/'atoms', '*')
        atmos.read_atom()
        departure = StudioDeparture.from_run(run)
        departure.atmos = atmos
        return departure