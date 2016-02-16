import numpy.random as ran
import numpy as np
import tasks as t
reload(t)

# number of rows and columns, number of ranks
n=30
ran.seed(4)

# generate the list of tasks
tasks=[]
k=0
for i in range(n):
    for j in range(n):
        tasks.append(t.task(i,j,n))
        k=k+1
ntask=k
im0=t.make_image(tasks)

# construct linked lists of neighbors for all tasks and count nbor types
def count_ranks(tasks):
    count=np.zeros(tasks[0].nrank,dtype=np.int32)
    for ij in range(np.size(tasks)):
        t.make_nbor_lists(tasks,ij)
        t.count_nbors(tasks[ij])
        rank=tasks[ij].rank
        count[rank]+=1
    return count
count=count_ranks(tasks)
vnb0=t.print_nbors(tasks)
print count

def dump1(task,old_rank,new_rank):
    t.count_nbors(task)
    print task.vnbors,task.ranks
    node=task.nbors.head
    i=task.i
    j=task.j
    while node:
        nbor=node.data
        print "{:4d}{:4d}{:4d}{:3b}{:2b}".format(
            nbor.i-i,nbor.j-j,nbor.rank,nbor.rank==old_rank,nbor.rank==new_rank)
        node=node.next
    t.count_nbors(task)
    print task.vnbors

def dump(task,old_rank,new_rank,gain,gain1):
    print "gains don't match, dumping:", gain, gain1, old_rank,new_rank
    print "initial state:"
    task.rank=old_rank
    dump1(task,old_rank,new_rank)
    print "final state:"
    task.rank=new_rank
    dump1(task,old_rank,new_rank)

do_print=False
def migrate(task,old_rank,new_rank,gain):
    gain1=task.vnbors
    task.rank=new_rank
    count[old_rank]-=1
    count[new_rank]+=1
    vnbors=t.count_nbors(task)
    gain1=2*(gain1-task.vnbors)
    #print gain,gain1
    if gain!=gain1:
        dump(task,old_rank,new_rank,gain,gain1)
    if do_print: print "{:4d}{:4d}{:4d}{:4d}{:4d}".format(task.ij,
        old_rank,new_rank,count[old_rank],count[new_rank])
    return vnbors

# Now take a shot at load balancing
def load_balance0(tasks,tolerance):
    total_gain=0
    average=tasks[0].n**2/tasks[0].nrank
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        t.count_nbors(task)
        rank=task.rank
        new_rank=task.hvnbors.argmin()
        gain=2*(task.hvnbors[rank]-task.hvnbors[new_rank])
        #print rank,new_rank,task.hvnbors[rank],task.hvnbors[new_rank],gain
        #print ranks[rank],ranks[max_rank],gain
        if gain>0 and count[rank]>average-tolerance and count[new_rank]<average+tolerance:
            migrate(task,rank,new_rank,gain)
            total_gain+=gain
    t.print_nbors(tasks)
    
# Now take a shot at load balancing
def load_balance1(tasks):
    total_gain=0
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        t.count_nbors(task)
        rank=task.rank
        ranks=task.ranks
        max_rank=ranks.argmax()
        gain=2*(ranks[max_rank]-ranks[rank])
        #print ranks[rank],ranks[max_rank],gain
        if gain>0 and count[rank]>count[max_rank]-(gain-1):
            migrate(task,rank,max_rank,gain)
            total_gain+=gain
        elif gain>2 and count[rank]>count[max_rank]-1:
            migrate(task,rank,max_rank,gain)
            total_gain+=gain
        elif gain>0 and count[rank]>count[max_rank]:
            migrate(task,rank,max_rank,gain)
            total_gain+=gain
    t.print_nbors(tasks)
    
# Now take a shot at load balancing
def load_balance2(tasks):
    total_gain=0
    # search for reduced virtual nbors, accepting slight load imbalance
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        t.count_nbors(task)
        ranks=task.ranks
        max_rank=ranks.argmax()
        gain=2*(ranks[max_rank]-ranks[task.rank])
        if gain>0 and count[task.rank]>count[max_rank]-(gain+2):
            migrate(task,task.rank,max_rank,gain)
            total_gain+=gain
    # improve load balance, while accepting slight increase in vnbors
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        t.count_nbors(task)
        ranks=task.ranks
        count0=task.n**2/task.nrank
        # search for best nbor rank to migrate to
        best=task
        best_gain=-9
        node=task.nbors.head
        while node:
            nbor=node.data
            if count[task.rank]>count0 and count[nbor.rank]<count0:
                gain=2*(ranks[nbor.rank]-ranks[task.rank])
                if gain>best_gain:
                    best=nbor
                    best_gain=gain
            node=node.next
        # accept if not too bad
        if best_gain>-3: 
            migrate(task,task.rank,best.rank,best_gain)
            total_gain+=best_gain
    #vnbors,nnbors=t.sum_nbors(tasks)
    #print "b:",vnbors,count
    
def iterate_load_balance(niter,tolerance):
    vnb0=t.sum_nbors(tasks)
    for iter in range(niter): 
        load_balance0(tasks,tolerance)
        vnb1=t.sum_nbors(tasks)
        if vnb1==vnb0: break
        vnb0=vnb1
    print count
    vnb1=t.print_nbors(tasks)

niter=20
iterate_load_balance(niter,1)
iterate_load_balance(niter,2)
iterate_load_balance(niter,3)
t.show_image(tasks)