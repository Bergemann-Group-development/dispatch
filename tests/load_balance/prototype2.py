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
print count
vnb0=t.print_nbors(tasks)
count=count_ranks(tasks)
print count

# Now take a shot at load balancing
def load_balance(tasks):
    for ij in range(np.size(tasks)):
        task=tasks[ij]
        rank=task.rank
        vnbors=tasks[ij].vnbors
        max_rank=task.ranks.argmax()
        if count[rank]>count[max_rank]:
            print "{:4d}{:4d}{:4d}{:4d}{:4d}".format(ij,
                rank,max_rank,count[rank],count[max_rank])
            task.rank=max_rank
            count[rank]-=1
            count[max_rank]+=1
            t.count_nbors(task)

for i in range(10): load_balance(tasks)
print count
vnb1=t.print_nbors(tasks)
print vnb0,vnb1