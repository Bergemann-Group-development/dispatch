from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .LineMask import LineMask
from .synth import run_synth
import ipyformkit as ifk
import tumult
from tumult.linelist import LineList
import copy

fields = {('rerun', 'use ga', 'use gv', 'use gq'): (None, True, True, True)}
rerun_form = ifk.Form(fields)
rerun_form.vbox.add_class('rerun-form')
rerun_button = rerun_form.widgets_dict['rerun'].wid

def init_rerun_form(studio):
    def rerun_on_click(btn):
        run = studio.browser.selected.object
        studio.plot_synth(rerun_synth, run)
        
    rerun_button.on_click(rerun_on_click)

def rerun_synth(studio, name, old_run):
  default = old_run.folder
  line_mask = old_run.line_mask
  llist = LineList(old_run.llist_html.get_dataframe())
  values = rerun_form.check_and_return_values()
  settings = {'composition_params':{
    'disable_ga': not values['use ga'],
    'disable_vdw': not values['use gv'],
    'disable_stark': not values['use gq'],
    }}
  
  if hasattr(old_run, 'fit_params'):
    if isinstance(old_run, tumult.CoG) and 'abnd' in old_run.fit_params:
      settings['composition_params']['abundances'] = f'{old_run.element_id}={old_run.fit_params["abnd"]:.2f}'

  run = run_synth(studio.synth_folder, name, settings, line_mask, default, llist, progress=studio.progress)
  run.llist_html.reset_html()
  if hasattr(old_run, 'fit_params'):
    run.fit_params = copy.deepcopy(old_run.fit_params)
  return run
  

