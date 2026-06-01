from bokeh.models import ColumnDataSource
from . import ServerManager
from .ServerManager import reset_xaxis, update_source_data
from .auxfuncs import *
from .synth_form import synth_form
import numpy as np
import ipywidgets as widgets
from .settings import saved_settings
from tumult.auxfuncs import get_weq
import ipyformkit as ifk

fields = {'method': ('trapz', 'gauss', 'voigt', 'asym_gauss'),
          ('EW [mÅ]', 'reduced'): (0.0, 0.0)}

tooltips = {'method': """trapz yields EW from trapzoidal integration only in the selected region.
The other methods fit a profile to the data in the selected region and integrate analytically."""}

formatters = {'EW [mÅ]': lambda x: f'{x:.3f}',
              'reduced': lambda x: f'{x:.3f}'}

disable = {'EW [mÅ]': True, 'reduced': True}
ew_form = ifk.Form(fields, disable=disable, tooltips=tooltips, formatters=formatters)
ew_method = ew_form.inputs_dict['method']

def init_ew_button(studio):
    new_ew_button = widgets.Button(description='Equivalent Width')
    new_ew_button.layout.display = 'none'
    
    @log_callback_errors
    def selection_callback(attr, old, new):
        nothing_selected = len(studio.fig_source.selected.indices)==0
        display = 'none' if nothing_selected else 'block'
        new_ew_button.layout.display = display
        
    @log_callback_errors
    def new_ew(btn):
        ews = studio.browser.get_objects_by_type(EquivalentWidth)
        color_cycle = ColorCycle(saved_settings['bokeh palette'])
        color = color_cycle.get_new_color([ew.color for ew in ews])
        new_ew = EquivalentWidth(studio.fig_source, studio.obs.source.data, color=color)
        label = f'{new_ew.x0:.1f}'
        studio.browser.add_object(new_ew, label=label)
        new_ew.draw(studio.fig)

    @log_callback_errors
    def on_method_change(change):
        ew = studio.current
        ew.get_equivalent_width()
    
    studio.fig_source.selected.on_change('indices', selection_callback)
    new_ew_button.on_click(new_ew)
    ew_method.observe(on_method_change, names='value')
    studio.controls.children += (new_ew_button,)

class EquivalentWidth():
    icon = 'chart-area'

    def __init__(self, source, data, color=None):
        islc = np.sort(np.array(source.selected.indices))
        self.xdata = source.data['x'][islc]
        self.ydata = source.data['y'][islc]
        self.xmin = self.xdata[0]
        self.xmax = self.xdata[-1]
        self.x0 = 0.5 * (self.xmin + self.xmax)
        self.color = color
        self.glyphs = []

        data = {'x': self.xdata, 'y2': np.ones_like(self.xdata)}
        self.source = ColumnDataSource(data)

    def chosen(self, studio):
        self.update_ew_fields()
        instructions = {'xmin': self.xmin,
                        'xmax': self.xmax,
                        'tabs': (ew_form,),
                        'titles': ('EW Settings',)}
        return instructions

    def draw(self, *args):
        @log_callback_errors
        def draw_varea():
            self.get_equivalent_width()

            for p in args:
                varea = p.varea('x', 'y1', 'y2', source=self.source, hatch_pattern='right_diagonal_line',
                                hatch_scale=10, hatch_alpha=0.5, fill_color=self.color, fill_alpha=0.1)
                line = p.line('x', 'y1', source=self.source, line_color=self.color)
                self.glyphs.append((p, varea))
                self.glyphs.append((p, line))

        ServerManager.server.next_tick_update(draw_varea)

    def get_equivalent_width(self):
        self.method = ew_form['method']
        self.weq, func = get_weq(self.xdata, self.ydata, lam0=self.x0, method=self.method, get_func=True)
        self.ewr = np.log10(self.weq * 1e3 / self.x0)
        self.update_ew_fields()
        
        @update_source_data(self.source)
        def update_data(data):
            data['y1'] = func(self.xdata)

        update_data()

    def update_ew_fields(self):
        ew_form['method'] = self.method
        ew_form['EW [mÅ]'] = self.weq
        ew_form['reduced'] = self.ewr

    def set_visbility(self, visible):
        def setting_visibility():
            for _, g in self.glyphs:
                g.visible = visible
        ServerManager.server.next_tick_update(setting_visibility)

    def delete(self):
        self.set_visbility(False)

    def hide(self):
        self.set_visbility(False)

    def show(self):
        self.set_visbility(True)

    def __str__(self):
        return self.label
