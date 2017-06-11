import matplotlib.pyplot as pl
import numpy as np

import dispatcher; reload(dispatcher)
from dispatcher import dispatcher

save=0
pl.figure(1,figsize=(10,5))
if not save: pl.show()

d=dispatcher()

node=d.tl.first
i=0
j=0
tasks=[]
while node and j<50:
    next=node.next
    task=d.tl.remove(node)
    task.update()
    tasks.append(task)
    d.add_by_time(task)
    node=next
    i=i+1
    ''' every nthreads time, make new image '''
    if np.mod(i,d.node.nthreads)==0:
        i=0
        if save:
            pl.close()  
            pl.figure(1,figsize=(10,5))
        pl.clf()
        pl.xlim(0,d.size[0])
        pl.ylim(0,d.size[1])
        d.show()
        ''' show the last nthreads tasks in mid color '''
        for task in tasks:
            task.show(color=0.5)
        tasks=[]
        pl.draw()
        pl.tight_layout()
        if save:
            file='anim/{:04d}.png'.format(j)
            pl.savefig(file)
            print file
        else:
            pl.pause(0.001)
        j=j+1