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
from matplotlib.mlab import find

# number of rows and columns, number of ranks
n=30
nrank=5
ran.seed(4)

rank=ran.random_integers(0,nrank-1,(n,n))
rank0=np.copy(rank)
count=np.zeros(nrank,dtype=np.int32)
ranks=np.zeros((n,n,nrank),dtype=np.int32)   # ranks[i,j,k] = number of [i,j] nbors with rank==k
tmp=np.zeros(nrank,dtype=np.int32)
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
            if i1!=0 or j1!=0: tmp[rank[i2,j2]]+=1
            if rank[i2,j2]!=rank[i,j]:
                nvn+=1
    return nvn

# Compute the number of virtual nbors for all tasks
def count_vnbors():
    count[:]=0
    for i0 in range(n):
        for j0 in range(n):
            tmp[:]=0
            count[rank[i0,j0]]+=1
            vnbors[i0,j0]=nvnb(i0,j0)
            ranks[i0,j0,:]=tmp[:]
    return np.sum(vnbors)
def print_vnbors():
    vnb=count_vnbors()
    print vnb, count

print_vnbors()

def load_balance(iter,tol):
    for i0 in range(n):
        for j0 in range(n):
            # the [i1,j1] index pair runs over all neighbors of [i0,j0] (including that point!)
            rank0=rank[i0,j0]
            # ranks[i0,,jo,k] is the number of nbors with rank=k
            nranks=ranks[i0,j0,:]
            choices=find(nranks!=0)                 # choices is list of ranks
            choices=choices[find(choices!=rank0)]   # exclude rank0
            nc=np.size(choices)                     # number of choices
            if nc==0: continue
            dvn=np.zeros(nc,dtype=np.int32)         # change in vnbors
            for i in range(nc):
                rank2=choices[i]
                if rank2==rank0: continue
                rank[i0,j0]=rank2
                # this double loop computes the change of the number of virtual 
                # neighbors in the 9 cells centered on [i0,j0]. 
                dn=0
                for i3 in range(-1,2):
                    i4=np.mod(i0+i3,n)
                    for j3 in range(-1,2):
                        j4=np.mod(j0+j3,n)
                        dn=dn+nvnb(i4,j4)-vnbors[i4,j4]
                # store, for this choice, the chance in vnb
                dvn[i]=dn
                rank[i0,j0]=rank0
            rank2=-1
            # if there is overload, and change that improves vnb
            if count[rank0]>count0:
                rank1=choices[dvn.argmin()]
                if dvn.min()<0:
                    if count[rank1]<count0:
                        rank2=rank1
                        print iter,i0,j0,"improving vnb and balance: change",rank0,rank2
            if rank2<0 and count[rank0]>count0-tol:
                rank1=choices[dvn.argmin()]
                if dvn.min()<0:
                    if count[rank1]<count0+tol:
                        rank2=rank1
                        print iter,i0,j0,"improving vnb",rank0,rank2
            # if no good case found, check that dvn increase with at most 1, and it is
            # possible to improve the balance
            if rank2<0 and dvn.min()<2:
                rank1=choices[dvn.argmin()]
                if count[rank0]>count0 and count[rank1]<count0:
                    rank2=rank1
                    print iter,i0,j0,"improving balance",rank0,rank2
            if rank2>-1:
                rank[i0,j0]=rank2
                count[rank0]-=1
                count[rank2]+=1
                for i3 in range(-1,2):
                    i4=np.mod(i0+i3,n)
                    for j3 in range(-1,2):
                        j4=np.mod(j0+j3,n)
                        vnbors[i4,j4]=nvnb(i4,j4)
                vnb=count_vnbors()
                print i0,j0,vnb,count
            #else:
                #print i0,j0,"no change"
    print_vnbors()
    #print xx0
    show_ranks(rank,iter)
    
def iterate_load_balance(niter,tolerance):
    vnb0=count_vnbors()
    for iter in range(1,niter+1):
        load_balance(iter,tolerance)
        vnb1=count_vnbors()
        if vnb1==vnb0: break
        vnb0=vnb1

niter=60
iterate_load_balance(niter,1)
iterate_load_balance(niter,2)
iterate_load_balance(niter,3)
iterate_load_balance(niter,4)

rank2=ran.random_integers(0,nrank-1,(n*2,n*2))
for i in range(2):
    for j in range(2):
        rank2[n*i:n*i+n,n*j:n*j+n]=rank

#pl.figure()
#show_ranks(rank2)
