from .auxfuncs import *
from .ServerManager import *
import ipyformkit as ifk
import ipywidgets as widgets
import numpy as np
import m3dis
import json


fields = {'synthesis folder': 'TumultStudio/tmp/',
          'atmospheres folder':'input_multi3d/atmos/',
          'bokeh palette': 'Light',
          'save':None}

with open(module_dir + 'settings.json', 'r') as f:
    saved_settings = json.load(f)
    
def save_settings(btn):
    new_settings = settings_form.check_and_return_values()
    saved_settings.update(new_settings)
    with open(module_dir + 'settings.json', 'w') as f:
        json.dump(saved_settings, f, indent=4)

settings_form = ifk.Form(fields)
settings_form.set_values(saved_settings)
save_settings_button = settings_form.widgets_dict['save'].wid
save_settings_button.on_click(save_settings)

title = widgets.Label('Global Settings')
close_button = widgets.Button(icon='times')
header = widgets.HBox([title, close_button])
header.add_class('header')
content = settings_form.vbox
settings_page = widgets.VBox([header, content])
settings_page.add_class('tum-settings')
settings_page.close_button = close_button
    
    