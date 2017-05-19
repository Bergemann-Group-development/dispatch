import numpy as np
import llist as ll
import task
reload(task)
import image
reload(image)

class thread:
    def __init__(self,id):
        self.id=id
        self.color=0.1*(1+id)

class dispatcher:
    size=[20,10]
    nthreads=6
    tl=ll.dllist()

    def add_by_time(self, task):
        if not self.tl.last:
            #print "first element",task.id
            self.tl.append(task)
        if task.t>=self.tl.last.value.t:
            #print "appending",task.id
            self.tl.append(task)
            return          
        if task.t<self.tl.first.value.t:
            #print "prepending",task.id
            self.tl.appendleft(task)
            return          
        node=self.tl.first
        while node:
            # print task.id,node.value.t, task.t
            if node.value.t>task.t:
                self.tl.insert(task,node)
                return
            node=node.next

    def __init__(self):
        id=0
        self.im=image.image(size=self.size)
        self.tl.clear()
        # create a task list
        for ix in range(self.size[0]):
            for iy in range(self.size[1]):
                llc=[ix,iy]
                size=[1,1]
                t=task.task(id,llc,size)
                if np.abs(iy-5.5)<1:
                    t.dt=t.dt/5.
                id=id+1
                self.add_by_time(t)
        # create a list of threads
        self.threads=[]
        for i in range(self.nthreads):
            self.threads.append(thread(i))
        
    def info(self):
        t=self.tl.first
        while t:
            t.value.info()
            t=t.next
