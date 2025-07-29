from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
from .LineMask import LineMask
from .synth_form import run_synth
import ipyformkit as ifk
import m3dis
from m3dis.linelist import LineList
import copy

fields = {('rerun', 'use ga', 'use gv', 'use gq'): (None, True, True, True)}
rerun_form = ifk.Form(fields)
rerun_form.vbox.add_class('rerun-form')

def rerun_synth(name, old_run, progress=None):
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
    if isinstance(old_run, m3dis.CoG) and 'abnd' in old_run.fit_params:
      settings['composition_params']['abundances'] = f'{old_run.element_id}={old_run.fit_params["abnd"]:.2f}'

  run = run_synth(name, settings, line_mask, default, llist, progress=progress)
  run.llist_html.reset_html()
  if hasattr(old_run, 'fit_params'):
    run.fit_params = copy.deepcopy(old_run.fit_params)
  return run
  

