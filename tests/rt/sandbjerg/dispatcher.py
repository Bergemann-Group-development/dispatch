''' Dispatcher class (was in task_list) '''
import numpy as np
import llist as ll

import task; reload(task)
import image; reload(image)
from task  import task
from node  import node
from image import image

class dispatcher:
    def add_by_time(self, task):
        ''' add a task, sorted by time '''
        node=self.tl.first
        if node:
            while node:
                if node.value.t>task.t:
                    self.tl.insert(task,node)
                    return
                node=node.next
        self.tl.append(task)
    def __init__(self):
        ''' initialize a dispatcher, with task list '''
        self.size=[20,10]
        self.im=image(size=self.size)
        self.node=node()
        self.tl=ll.dllist()
        self.tl.clear()
        id=0
        for ix in range(self.size[0]):
            for iy in range(self.size[1]):
                t=task(id,llc=[ix,iy],size=[1,1])
                if np.abs(iy-5.5)<1:
                    t.dt=t.dt/5.
                id=id+1
                self.add_by_time(t)
    def show(self,csz=1.0):
        ''' show all tasks with colored time '''
        item=self.tl.first
        while item:
            task=item.value
            task.show(color=1.0-np.mod(task.t,csz)/csz)  
            item=item.next
