import tumult
from .ServerManager import update_data_source, clear_data
from .nlte_form import nlte_form
from .auxfuncs import *
from .synth import StudioSynth
import ipyformkit as ifk
import numpy as np

mhd_vars = ('ltau', 'temp', 'rho', 'ne', 'chi', 'vx', 'vy', 'vz')

fields = {
    'source': ('Atmosphere',),
    ('x-axis', 'y-axis'): (mhd_vars, mhd_vars),}

atmos_plot_form = ifk.Form(fields)
source_wid = atmos_plot_form.widgets_dict['source'].wid
xquantity_wid = atmos_plot_form.widgets_dict['x-axis'].wid
yquantity_wid = atmos_plot_form.widgets_dict['y-axis'].wid

def init_atmos_plot_form(studio):
    @log_callback_errors
    def on_source_change(change):
        atmos = studio.browser.selected.object
        val = source_wid.value
        if val == 'Atmosphere':
            xquantity_wid.options = mhd_vars
            yquantity_wid.options = mhd_vars
        else:
            _, dep = atmos.get_depart(atom=val, ilvl=None)
            nlvl = dep.shape[-1]
            xquantity_wid.options = mhd_vars
            yquantity_wid.options = np.arange(nlvl)

    @log_callback_errors
    def on_quantity_change(change):
        atmos = studio.browser.selected.object
        atmos.draw_atmos_data(studio)
        
    xquantity_wid.value = 'ltau'
    yquantity_wid.value = 'temp'
    source_wid.observe(on_source_change, names='value')
    xquantity_wid.observe(on_quantity_change, names='value')
    yquantity_wid.observe(on_quantity_change, names='value')

class StudioAtmos(StudioSynth):
    def add_plot_source(self, source):
        source_wid.options += (source,)

    def draw_atmos_data(self, studio):
        val = source_wid.value
        @update_data_source(studio.atmos_source)
        def update_atmos_data(data):
            if val == 'Atmosphere':
                data['x'] = self.__dict__[xquantity_wid.value].ravel()
                data['y'] = self.__dict__[yquantity_wid.value].ravel()
                studio.atmos_fig.xaxis.axis_label = xquantity_wid.value
                studio.atmos_fig.yaxis.axis_label = yquantity_wid.value
            else:
                ilvl = yquantity_wid.value
                data['x'] = self.__dict__[xquantity_wid.value]
                studio.atmos_fig.xaxis.axis_label = xquantity_wid.value
                _, dep = self.get_depart(atom=val, ilvl=ilvl)
                data['y'] = dep[..., ::-1]
                studio.atmos_fig.yaxis.axis_label = f'b_{ilvl}'
        update_atmos_data()
    
    def chosen(self, studio):
        self.draw_atmos_data(studio)
        studio.right_sidebar.children = (nlte_form.vbox, atmos_plot_form.vbox)
        studio.right_sidebar.titles = ['NLTE Settings', 'Plot Settings']