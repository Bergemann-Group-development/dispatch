'''
   1) Make a square (or cubic) array, with NxN tasks
   2) Assign a rank randomly to each position
   3) For each "task", compute the number of neighbors w different rank
   4) For each neighbor of the task, compute the increase or decrease in
      the number of virtual tasks, if the rank of the local task was
      changed to that of the neighbor
'''
import numpy as np
import numpy.random as ran
import matplotlib.pyplot as pl

# number of rows and columns, number of ranks
n=30
nrank=5
ran.seed(3)

rank=ran.random_integers(0,nrank-1,(n,n))
rank0=np.copy(rank)
count=np.zeros(nrank,dtype=np.int32)
count0=n*n/nrank

pl.ion()
def show_ranks(r,iter):
    pl.imshow(r,interpolation='nearest')
    pl.title(iter)
    pl.draw()
    pl.pause(0.1)
    
pl.figure(1,figsize=(8,8))
show_ranks(rank0,0)

#nbors =np.zeros((n,n),dtype=np.int32)
vnbors=np.zeros((n,n),dtype=np.int32)

# Count the number of virtual neighbors for one task
def nvnb(i,j):
    nvn=0
    for i1 in range(-1,2):
        i2=np.mod(i+i1,n)
        for j1 in range(-1,2):
            j2=np.mod(j+j1,n)
            if rank[i2,j2]!=rank[i,j]:
                nvn+=1
    return nvn

# Compute the number of virtual nbors for all tasks
def count_vnbors():
    count[:]=0
    for i0 in range(n):
        for j0 in range(n):
            count[rank[i0,j0]]+=1
            vnbors[i0,j0]=nvnb(i0,j0)
    print np.sum(vnbors), count
    #print count,count0

count_vnbors()

def load_balance(iter,tol):
    for i0 in range(n):
        for j0 in range(n):
            for i1 in range(-1,2):
                i2=np.mod(i0+i1,n)
                for j1 in range(-1,2):
                    j2=np.mod(j0+j1,n)
                    rank0=rank[i0,j0]
                    rank2=rank[i2,j2]
                    # if rank2 is different from rank0, then a change of
                    # rank0 to rank2 means the number of vnbors[i2,j2] is 
                    # reduced by one, while vnbors[i0,j0] must be recomputed
                    if rank2!=rank0 and count[rank2]<(count0+tol) and count[rank0]>(count0-tol):
                        rank[i0,j0]=rank2
                        dn=0
                        for i3 in range(-1,2):
                            i4=np.mod(i0+i3,n)
                            for j3 in range(-1,2):
                                j4=np.mod(j0+j3,n)
                                dn=dn+nvnb(i4,j4)-vnbors[i4,j4]
                        # if the number of vnbors is increased by a migration, cancel
                        if dn>0:
                            rank[i0,j0]=rank0
                        # else update rank and vnb counts
                        else:
                            count[rank2]+=1
                            count[rank0]-=1
                            for i3 in range(-1,2):
                                i4=np.mod(i0+i3,n)
                                for j3 in range(-1,2):
                                    j4=np.mod(j0+j3,n)
                                    vnbors[i4,j4]=nvnb(i4,j4)
                            #print i0,j0,vnbors[i0,j0],nvnb(i0,j0),vnbors[i2,j2],nvnb(i2,j2),dvnbs
    count_vnbors()
    show_ranks(rank,iter)
    
def iterate_load_balance(niter,tolerance):
    for iter in range(1,niter+1):
        load_balance(iter,tolerance)

niter=20
iterate_load_balance(niter,1)
iterate_load_balance(niter,2)
iterate_load_balance(niter,3)

rank2=ran.random_integers(0,nrank-1,(n*2,n*2))
for i in range(2):
    for j in range(2):
        rank2[n*i:n*i+n,n*j:n*j+n]=rank

#pl.figure()
#show_ranks(rank2)
