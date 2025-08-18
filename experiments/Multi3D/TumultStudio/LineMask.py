from bokeh.models import *
from . import ServerManager
from .ServerManager import reset_xaxis
from .auxfuncs import *
from .synth_form import synth_form
import numpy as np
import ipywidgets as widgets
from .settings import saved_settings

def init_line_mask_button(studio):
    new_line_mask_button = widgets.Button(description='New Line Mask')
    new_line_mask_button.layout.display = 'none'
    
    @log_callback_errors
    def selection_callback(attr, old, new):
        nothing_selected = len(studio.fig_source.selected.indices)==0
        display = 'none' if nothing_selected else 'block'
        new_line_mask_button.layout.display = display
        
    @log_callback_errors
    def new_line_mask(btn):
        line_masks = studio.browser.get_objects_by_type(LineMask)
        color_cycle = ColorCycle(saved_settings['bokeh palette'])
        color = color_cycle.get_new_color([lm.color for lm in line_masks])
        new_mask = LineMask(studio.fig_source, studio.obs.source.data, color=color)
        new_mask.label = f'{new_mask.xmin:.1f}-{new_mask.xmax:.1f}'
        studio.browser.add_object(new_mask)
        new_mask.draw_boxes(studio.nav, studio.fig)
    
    studio.fig_source.selected.on_change('indices', selection_callback)
    new_line_mask_button.on_click(new_line_mask)
    studio.controls.children += (new_line_mask_button,)

class LineMask():
    def __init__(self, source, data, color=None):
        xx = source.data['x']
        islc = np.sort(np.array(source.selected.indices))
        splits = np.where(np.diff(islc)>1)[0]+1
        regions = np.split(islc, splits)
        self.xmin = xx[islc[0]]
        self.xmax = xx[islc[-1]]

        self.regions = []
        self.xdata = []
        self.ydata = []
        for reg in regions:
            xmin = xx[reg[0]]
            xmax = xx[reg[-1]]
            self.regions.append((xmin, xmax))
            slc = slicesorted(data['x'], xmin, xmax)
            self.xdata = np.concatenate([self.xdata, data['x'][slc]])
            self.ydata = np.concatenate([self.ydata, data['y'][slc]])

        self.boxes = []
        self.cntm = np.ones_like(self.xdata)
        self.color = color

    def chosen(self, studio):
        instructions = {'xmin': self.xmin,
                        'xmax': self.xmax,
                        'tabs': (synth_form.vbox,),
                        'titles': ('Synthesis Settings',)}
        return instructions

    def normalise(self, spline):
        # Line Mask data is always normalised
        # If renormalised, the old normalisation is undone
        self.ydata *= self.cntm
        self.cntm = spline(self.xdata)
        self.ydata /= self.cntm

    def get_wavelength_grid(self, daa=0.01):
        xdata = []
        for xmin, xmax in self.regions:
            n = np.ceil((xmax-xmin)/daa).astype(int) + 1 
            xdata = np.concatenate([xdata, np.linspace(xmin, xmax, n)])
        return xdata

    def draw_boxes(self, *args):
        for p in args:
            for xmin, xmax in self.regions:
                box = BoxAnnotation(left=xmin, right=xmax, line_color=None,
                                    fill_color=self.color, fill_alpha=0.1)
                self.boxes.append((p, box))
            
        def update():
            [p.add_layout(box) for p, box in self.boxes]
            
        ServerManager.server.next_tick_update(update)

    def delete(self):
        def update():
            for _, box in self.boxes:
                    box.visible = False
        ServerManager.server.next_tick_update(update)

    def hide(self):
        def update():
            for _, box in self.boxes:
                    box.visible = False
        ServerManager.server.next_tick_update(update)

    def show(self):
        def update():
            for _, box in self.boxes:
                    box.visible = True
        ServerManager.server.next_tick_update(update)

    def __str__(self):
        return self.label
