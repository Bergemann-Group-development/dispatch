from bokeh.models import *
from .ServerManager import *
import numpy as np

class LineMask():
    def __init__(self, source, label='untitled'):
        self.source = source
        self.islc = np.sort(np.array(source.selected.indices))
        self.splits = np.where(np.diff(self.islc)>1)[0]+1
        self.regions = np.split(self.islc, self.splits)
        self.label = label
        self.boxes = []

    def draw_boxes(self, *args, **kwargs):
        for p in args:
            for reg in self.regions:
                xmin = self.source.data['x'][reg[0]]
                xmax = self.source.data['x'][reg[-1]]
                box = BoxAnnotation(left=xmin, right=xmax, **kwargs)
                self.boxes.append((p, box))
            
        def update():
            [p.add_layout(box) for p, box in self.boxes]
            
        server_manager.next_tick_update(update)

    def jump_to_mask(self, p):
        xmin = self.source.data['x'][self.islc[0]]
        xmax = self.source.data['x'][self.islc[-1]]

        def update():
            p.x_range.start = xmin
            p.x_range.end = xmax
        
        server_manager.next_tick_update(update)

    def delete(self):
        def update():
            for _, box in self.boxes:
                    box.visible = False
        server_manager.next_tick_update(update)

    def __str__(self):
        return self.label
