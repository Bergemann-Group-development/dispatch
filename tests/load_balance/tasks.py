import lists as l
import numpy as np
import numpy.random as ran
import matplotlib.pyplot as pl

# class definition for tasks
class task:
    nrank=5
    def __init__(self,i,j,n):
        self.i=i
        self.j=j
        self.ij=j+n*i
        self.n=n
        self.rank=ran.random_integers(0,self.nrank-1)
        self.nnbors=0
        self.vnbors=0

# translate integer tuple to linear tasks index
def index(i,j,n):
    return j+i*n

# construct linked lists of neighbors for a task
def make_nbor_lists(tasks,ij):
    task=tasks[ij]
    task.nbors=l.linked_list()
    for i1 in range(-1,2):
        i2=np.mod(task.i+i1,task.n)
        for j1 in range(-1,2):
            if i1!=0 or j1!=0:
                j2=np.mod(task.j+j1,task.n)
                ij2=index(i2,j2,task.n)
                task.nbors.insert_node(tasks[ij2])

# count neighbors of a task
do_print=False
def count_nbors(task):
    task.nnbors=0
    task.vnbors=0
    node=task.nbors.head
    task.ranks=np.zeros(task.nrank,dtype=np.int32)
    while node:
        nbor=node.data
        task.ranks[nbor.rank]+=1
        if do_print: print "{0:3d}{1:3d}{2:3d}{3:3d}{4:3d}".format(
          task.j+task.i*task.n,
          np.mod(task.i-nbor.i+1,task.n)-1,
          np.mod(task.j-nbor.j+1,task.n)-1,
          task.rank,nbor.rank)
        if nbor.rank==task.rank:
            task.nnbors+=1
        else:
            task.vnbors+=1
        node=node.next

# print the number of vnbors
def sum_nbors(tasks):
    vnbors=0
    nnbors=0
    for ij in range(np.size(tasks)):
        count_nbors(tasks[ij])
        vnbors+=tasks[ij].vnbors
        nnbors+=tasks[ij].nnbors
    return (vnbors,nnbors)
    
# print the number of vnbors
def print_nbors(tasks):
    vnbors,nnbors=sum_nbors(tasks)
    print vnbors
    return vnbors

def make_image(tasks):
    n=tasks[0].n
    im=np.zeros((n,n),dtype=np.int32)
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        i=task.i
        j=task.j
        im[i,j]=task.rank
    return im
    
def show_image(tasks):
    im=make_image(tasks)
    pl.figure(1,figsize=(6,6))
    pl.imshow(im,interpolation='nearest')
    pl.colorbar()
    pl.show()