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
            changes=0
            dnmin=0
            rank0=rank[i0,j0]
            # ranks[i0,,jo,k] is the number of nbors with rank=k
            nranks=ranks[i0,j0,:]
            choices=find(nranks!=0)
            for rank2 in choices:
                if rank2==rank0: continue
                # if rank2 is different from rank0, 
                if count[rank2]<(count0+tol) and count[rank0]>(count0-tol):
                    rank[i0,j0]=rank2
                    count[rank2]+=1
                    count[rank0]-=1
                    dn=0
                    # this double loop computes the change of the number of virtual 
                    # neighbors in the 9 cells centered on [i0,j0]. 
                    for i3 in range(-1,2):
                        i4=np.mod(i0+i3,n)
                        for j3 in range(-1,2):  # [i4,j4] is the neighbor-neighbor
                            j4=np.mod(j0+j3,n)
                            dn=dn+nvnb(i4,j4)-vnbors[i4,j4]
                    # If the number of vnbors is increased by a migration, cancel!
                    if dn>0:
                        rank[i0,j0]=rank0
                        count[rank2]-=1
                        count[rank0]+=1
                    # else update rank and vnb counts
                    else:
                        #dnlast=dn
                        if dn<dnmin: 
                            changes+=1
                            dnmin=dn
                            rank1=rank2
                        # the [i0,j0] cell is changed from rank0 to rank2
                        #count[rank2]+=1
                        #count[rank0]-=1
                        #ranks[i0,i0,rank0]-=1
                        #ranks[i0,i0,rank2]+=1
                        # this loops over all nbors of [i0,j0], updating the count of
                        # virtual nbors of each nbor
                        #for i3 in range(-1,2):
                        #    i4=np.mod(i0+i3,n)
                        #    for j3 in range(-1,2):
                        #        j4=np.mod(j0+j3,n)
                        #        vnbors[i4,j4]=nvnb(i4,j4)
                        #xx0+=dn
                        #s="migrate:{:4d}{:4d}{:4d}{:4d}"
                        #print ranks[i0,j0,:],rank0,rank2
                        #print s.format(i0,j0,dn,dn1) #,dn2
                        #print i0,j0,vnbors[i0,j0],nvnb(i0,j0),vnbors[i2,j2],nvnb(i2,j2),dvnbs
            if changes>0: # and dnlast!=dnmin: 
                #print "changes={} last choice={} dn={} best={}".format(changes,rank2,dnlast,dnmin)
                rank[i0,j0]=rank1
                count[rank2]-=1
                count[rank1]+=1
                for i3 in range(-1,2):
                    i4=np.mod(i0+i3,n)
                    for j3 in range(-1,2):
                        j4=np.mod(j0+j3,n)
                        vnbors[i4,j4]=nvnb(i4,j4)
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

niter=50
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
