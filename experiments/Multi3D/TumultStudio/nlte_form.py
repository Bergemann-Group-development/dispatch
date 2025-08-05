from .ServerManager import update_data_source, clear_data
from .settings import saved_settings
from .synth import run_synth
from .auxfuncs import *
import ipyformkit as ifk
import numpy as np
import os
import logging
import pathlib

atom_files = [f for f in os.listdir(saved_settings['atom folder']) if not f.startswith('.')]

fields = {'atom file': tuple(atom_files),
          ('abundance', 'relative'): (0.0, True),
          'Compute Departures': None,}

nlte_form = ifk.Form(fields)
compute_departures_btn = nlte_form.widgets_dict['Compute Departures'].wid
atom_file_wid = nlte_form.widgets_dict['atom file'].wid
departure_select = widgets.Select(options=['no departure coefficients yet!'], description='')
departure_select.layout.width = '100%'
nlte_form.vbox.children = (departure_select,) + nlte_form.vbox.children

def init_nlte_form(studio):
    @log_callback_errors
    def on_compute_click(btn):
        launch_nlte(studio)
        
    compute_departures_btn.on_click(on_compute_click)

  
def get_nlte_settings(atmos):
    values = nlte_form.check_and_return_values()
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
    

def launch_nlte(studio):
    logging.info('Launching NLTE computation...')
    atmos = studio.browser.selected.object
    default = module_dir / 'default_nlte.nml'
    settings, name = get_nlte_settings(atmos)
    if settings:
        studio.progress.description = 'Running NLTE calculation...'
        run = run_synth(studio.synth_folder, name, settings, default=default, progress=studio.progress)
        atomid = run.atoms[0].id
        link_files(run.sfolder, atmos.sfolder, 'pop_*')
        link_files(run.folder/'atoms', atmos.folder/'atoms', '*')
        if departure_select.options[0].startswith('no dep'):
            departure_select.options = []
        if atomid not in departure_select.options:
            departure_select.options += (atomid,)
        atmos.add_plot_source(atomid)
        studio.progress.value = studio.progress.max
        studio.progress.description = 'Done!'
        return run
    else:
        return None