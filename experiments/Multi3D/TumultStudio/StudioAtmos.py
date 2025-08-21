import tumult
from .ServerManager import update_source_data, clear_data
from .departure_form import departure_form
from .auxfuncs import *
from .widgets import HTMLDataFrame
from .synth import StudioSynth
import ipyformkit as ifk
import numpy as np
import pandas as pd

mhd_vars = ('ltau', 'temp', 'rho', 'ne', 'chi', 'vx', 'vy', 'vz')

fields = {
    ('x-axis', 'y-axis'): (mhd_vars, mhd_vars),}

atmos_plot_form = ifk.Form(fields)
xaxis_wid = atmos_plot_form.widgets_dict['x-axis'].wid
yaxis_wid = atmos_plot_form.widgets_dict['y-axis'].wid

def init_atmos_plot_form(studio):
    @log_callback_errors
    def on_quantity_change(change):
        yaxis_wid.options = [o for o in mhd_vars if o != xaxis_wid.value]
        atmos = studio.current
        atmos.draw_atmos_data(studio)
        
    xaxis_wid.value = 'ltau'
    yaxis_wid.value = 'temp'
    xaxis_wid.observe(on_quantity_change, names='value')
    yaxis_wid.observe(on_quantity_change, names='value')

class StudioAtmos(StudioSynth):
    icon = 'sun'

    def _init_source(self):
        pass

    def draw_atmos_data(self, studio):
        clear_data(studio.atmos_interp_source)
        xval = xaxis_wid.value
        yval = yaxis_wid.value

        xscat = self.__dict__[xval]
        yscat = self.__dict__[yval]

        @update_source_data(studio.atmos_scatter_source)
        def update_atmos_scatter(data):
            studio.atmos_scatter_source.selected.indices = []
            data['x'] = np.ravel(xscat)
            data['y'] = np.ravel(yscat)
            studio.atmos_fig.xaxis.axis_label = xval
            studio.atmos_fig.yaxis.axis_label = yval
        update_atmos_scatter()

    def chosen(self, studio):
        self.draw_atmos_data(studio)
        instructions = {'tabs':(departure_form.vbox, atmos_plot_form.vbox, self.info),
                        'titles':('NLTE Settings', 'Plot Settings', 'Info'),}
        return instructions

    @property
    def info(self):
        data = {'File': self.nml['atmos_params']['atmos_file'].strip(),
                'Format': self.nml['atmos_params']['atmos_format'].strip(),
                'FeH': self.nml['atmos_params']['feh'],
                'vmic': self.nml['atmos_params']['vmic'],}
        if 'nx' in self.nml['atmos_params']:
            data['nx'] = self.nml['atmos_params']['nx']
            data['ny'] = self.nml['atmos_params']['ny']
        df = pd.DataFrame(list(data.items()), columns=['key', 'value'])
        return HTMLDataFrame(df)