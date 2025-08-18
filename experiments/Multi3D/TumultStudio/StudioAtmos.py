import tumult
from .ServerManager import update_source_data, clear_data
from .departure_form import departure_form, set_departure_form
from .auxfuncs import *
from .widgets import HTMLDataFrame
from .synth import StudioSynth
import ipyformkit as ifk
import numpy as np
import pandas as pd

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

def set_atmos_plot_form(atmos):
    source_wid.options = ('Atmosphere',)
    source_wid.value = source_wid.options[0]

class StudioAtmos(StudioSynth):
    #=================================================================
    def _init_source(self):
        pass

    def add_plot_option(self, option):
        source_wid.options += (option,)

    def draw_atmos_data(self, studio):
        val = source_wid.value
        xval = xquantity_wid.value
        yval = yquantity_wid.value
        xlabel = xval

        xscat = self.__dict__[xval]
        if val == 'Atmosphere':
            yscat = self.__dict__[yval]
            ylabel = yval
        else:
            _, dep = self.get_depart(atom=val, ilvl=yval)
            yscat = dep[..., ::-1]
            ylabel = f'b_{yval}'

        @update_source_data(studio.atmos_scatter_source)
        def update_atmos_scatter(data):
            data['x'] = np.ravel(xscat)
            data['y'] = np.ravel(yscat)
            studio.atmos_fig.xaxis.axis_label = xlabel
            studio.atmos_fig.yaxis.axis_label = ylabel
        update_atmos_scatter()

        if xval == 'ltau':
            xline, yline = self.get_ltau_mean(yscat)
        else:
            xline = []
            yline = []

        @update_source_data(studio.atmos_line_source)
        def update_atmos_line(data):
            data['x'] = xline
            data['y'] = yline
        update_atmos_line()

    
    def chosen(self, studio):
        set_atmos_plot_form(self)
        set_departure_form(self)
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