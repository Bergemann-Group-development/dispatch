import matplotlib.pyplot as pl
import numpy as np
import task_list
reload(task_list)
from task_list import *

import image
reload(image)
from image import *

pl.figure(1)
d=dispatcher()

#d.info()
#d.im.rectangle(color=0.7)

node=d.tl.first
i=0
while node:
    next=node.next
    task=d.tl.remove(node)
    print task.id,task.t
    i=np.int(np.random.uniform(d.nthreads))
    task.update()
    task.show(d.im,color=d.threads[i].color)
    d.add_by_time(task)
    node=next
    i=i+1
    if i>=d.nthreads:
        i=0
        pl.draw()
        pl.pause(2)
        d.im.reset()
