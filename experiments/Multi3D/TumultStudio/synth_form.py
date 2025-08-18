from .auxfuncs import *
from .ServerManager import *
from .StudioAtmos import StudioAtmos
from .synth import run_synth
from .settings import saved_settings
import ipyformkit as ifk
import numpy as np

fields = {'atmosphere':('no atmospheres imported yet',),
          ('abundances', 'relative'):('O=8.70, Fe=7.50', False),
          ('exclude elements','only elements'):('O, Fe', 'O, Fe'),
          'launch':None}

tooltips = {
    'relative':'If true please specify abundances as [X/Fe], otherwise please specify absolute abundances X/H',
    'exclude elements':'Elements to exclude from the linelist.',
    'only elements':'All other elements will be removed from the linelist. If empty, all elements will be included.',
    }

disable = {
    'exclude elements': lambda d: d['only elements'] != '',
    }

check = {'atmosphere': lambda d: isinstance(d['atmosphere'], StudioAtmos),
         }

synth_form = ifk.Form(fields, tooltips=tooltips, disable=disable, check=check)
atmos_wid = synth_form.widgets_dict['atmosphere'].wid
launch_button = synth_form.widgets_dict['launch'].wid


def init_synth_form(studio):
    def launch_on_click(btn):
        studio.plot_synth(launch_synth)

    launch_button.on_click(launch_on_click)

    def update_atmos_dropdown(change):
        atmospheres = studio.browser.get_objects_by_type(StudioAtmos)
        if atmospheres:
            atmos_wid.options = tuple(atmospheres)
        else:
            atmos_wid.options = ('no atmospheres imported yet',)

    studio.browser.observe(update_atmos_dropdown, names='objects')
    update_atmos_dropdown(None)


def get_synth_settings():
    values = synth_form.check_and_return_values()
    if values:
        settings = {'atmos_params':{'atmos_file':values['atmosphere'].folder,
                                    'atmos_format':'save',
                                },
                    'composition_params':{'abundances':values['abundances'],
                                        'relative_abundances':values['relative'],
                                        },
                    'linelist_params':{'exclude_elements':values['exclude elements'],
                                    'only_elements':values['only elements']
                                    }}
        return settings
    else:
        return None


def launch_synth(studio, name, **kwargs):
    logging.info('Launching synth computation...')
    line_mask = studio.browser.selected.object
    atmos = atmos_wid.value 
    settings = get_synth_settings()
    if not settings:
        return None
    
    run = run_synth(studio, name, settings, line_mask, **kwargs)
    if run: run.atmos = atmos
    return run