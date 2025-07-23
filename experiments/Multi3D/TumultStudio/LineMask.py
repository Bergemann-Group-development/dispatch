from bokeh.models import *
from .ServerManager import *
from .auxfuncs import *
import numpy as np

class LineMask():
    def __init__(self, source, data, label='untitled'):
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

        self.label = label
        self.boxes = []

    def draw_boxes(self, *args, **kwargs):
        for p in args:
            for xmin, xmax in self.regions:
                box = BoxAnnotation(left=xmin, right=xmax, **kwargs)
                self.boxes.append((p, box))
            
        def update():
            [p.add_layout(box) for p, box in self.boxes]
            
        server_manager.next_tick_update(update)

    def jump_to_mask(self, p):
        def update():
            p.x_range.start = self.xmin
            p.x_range.end = self.xmax
        
        server_manager.next_tick_update(update)

    def delete(self):
        def update():
            for _, box in self.boxes:
                    box.visible = False
        server_manager.next_tick_update(update)

    def __str__(self):
        return self.label
