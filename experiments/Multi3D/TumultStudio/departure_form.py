from .ServerManager import update_source_data, clear_data
from .settings import saved_settings
from .synth import run_synth
from .auxfuncs import *
from . import ServerManager
from .threading import execute_in_thread
import ipyformkit as ifk
import numpy as np
import os
import logging
import pathlib

atom_files = [f for f in os.listdir(saved_settings['atom folder']) if not f.startswith('.')]

fields = {'atom file': tuple(atom_files),
          ('abundance', 'relative'): (0.0, True),
          'Compute Departures': None,
          'Cancel Departures': None}

departure_form = ifk.Form(fields)
compute_departures_btn = departure_form.widgets_dict['Compute Departures'].wid
cancel_departures_btn = departure_form.widgets_dict['Cancel Departures'].wid
cancel_departures_btn.layout.display = 'none'  # Hide cancel button initially
atom_file_wid = departure_form.widgets_dict['atom file'].wid
departure_select = widgets.Select(options=['no departure coefficients yet!'], description='')
departure_select.layout.width = '100%'
departure_form.vbox.children = (departure_select,) + departure_form.vbox.children

def init_departure_form(studio):
    @log_callback_errors
    def on_compute_click(btn):
        launch_departure(studio)
        
    compute_departures_btn.on_click(on_compute_click)

def set_departure_form(atmos):
    atomids = [atom.id for atom in atmos.atoms]
    if len(atomids) == 0:
        departure_select.options = ['no departure coefficients yet!']
    else:
        departure_select.options = atomids
        for atomid in atomids:
            atmos.add_plot_option(atomid)
  
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
    

def launch_departure(studio):
    container = studio.browser.selected
    atmos = container.object
    default = module_dir / 'default_nlte.nml'
    settings, name = get_departure_settings(atmos)
    if settings is None: return

    @execute_in_thread(cancel=cancel_departures_btn, spinner=container.label)
    def async_callback(execution):
        run = run_synth(studio, name, settings, default=default, runner=execution.run_command)
        if run:
            atomid = run.atoms[0].id
            link_files(run.sfolder, atmos.sfolder, 'pop_*')
            link_files(run.folder/'atoms', atmos.folder/'atoms', '*')
            atmos.read_atom()
            set_departure_form(atmos)

        studio.progress.value = studio.progress.max
        studio.progress.description = 'Done!'
        cancel_departures_btn.layout.display = 'none'
        cancel_departures_btn._click_handlers.callbacks.clear()
    
    studio.progress.description = 'Running departure calculation...'
    cancel_departures_btn.layout.display = 'block'  # Show cancel button
    async_callback()