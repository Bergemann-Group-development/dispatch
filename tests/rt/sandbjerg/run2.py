import matplotlib.pyplot as pl
import numpy as np
import image

from dispatcher import *

pl.figure(1,figsize=(2,1))
pl.clf()
#pl.show()

d=dispatcher()

node=d.tl.first
i=0
j=0
tasks=[]
while node and j<500:
    next=node.next
    task=d.tl.remove(node)
    #print task.id,task.t
    i=np.int(np.random.uniform(d.node.nthreads))
    task.update()
    tasks.append(task)
    d.add_by_time(task)
    node=next
    i=i+1
    if np.mod(i,d.node.nthreads)==0:
        pl.clf()
        tmp=d.tl.first
        csz=1.0
        while tmp:
            task=tmp.value
            task.show(color=1.0-np.mod(task.t,csz)/csz)  
            tmp=tmp.next
        for task in tasks:
            task.show(color=0.5)
        tasks=[]
        pl.xlim(0,d.size[0])
        pl.ylim(0,d.size[1])
        pl.draw()
        pl.close()
        file='anim/{:04d}.png'.format(j)
        pl.savefig(file)
        print file
        j=j+1
