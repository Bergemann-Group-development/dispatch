import tumult
from .ServerManager import update_source_data, clear_data
from .auxfuncs import *
from .widgets import HTMLDataFrame
from .synth import StudioSynth
import ipyformkit as ifk
import numpy as np
import pandas as pd

mhd_vars = ('ltau', 'temp', 'rho', 'ne', 'chi', 'vx', 'vy', 'vz')

fields = {
    ('x-axis', 'y-axis'): (mhd_vars, (0,)),}

departure_plot_form = ifk.Form(fields)
xaxis_wid = departure_plot_form.inputs_dict['x-axis']
yaxis_wid = departure_plot_form.inputs_dict['y-axis']

def init_departure_plot_form(studio):
    xaxis_wid.options = mhd_vars

    @log_callback_errors
    def on_quantity_change(change):
        departure = studio.current
        departure.draw_departure_data(studio)
        
    xaxis_wid.value = 'ltau'
    yaxis_wid.value = 0
    xaxis_wid.observe(on_quantity_change, names='value')
    yaxis_wid.observe(on_quantity_change, names='value')

class StudioDeparture(StudioSynth):
    icon = 'atom'

    def _init_source(self):
        pass
    
    def draw_departure_data(self, studio):
        clear_data(studio.atmos_interp_source)
        xval = xaxis_wid.value
        yval = yaxis_wid.value
        xlabel = xval
        ylabel = f'b_{yval}'

        xscat = self.__dict__[xval]
        _, dep = self.get_depart(atom=0, ilvl=yval)
        yscat = dep[..., ::-1]

        @update_source_data(studio.atmos_scatter_source)
        def update_atmos_scatter(data):
            studio.atmos_scatter_source.selected.indices = []
            data['x'] = np.ravel(xscat)
            data['y'] = np.ravel(yscat)
            studio.atmos_fig.xaxis.axis_label = xlabel
            studio.atmos_fig.yaxis.axis_label = ylabel
        update_atmos_scatter()

    def chosen(self, studio):
        yaxis_wid.options = np.arange(self.atoms[0].nk)
        self.draw_departure_data(studio)
        instructions = {'tabs':(departure_plot_form,),
                        'titles':('Plot Settings',),}
        return instructions