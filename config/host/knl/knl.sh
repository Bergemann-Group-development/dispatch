#!/bin/bash

# If not set, MPI will do it automagically to the right value (using all HW threads)
#export OMP_NUM_THREADS=16

# Decide how many resources each OMP thread should have; a core or a thread?
#export KMP_AFFINITY=verbose,scatter,granularity=thread
#export KMP_AFFINITY=verbose,scatter,granularity=core
#export KMP_AFFINITY=scatter,granularity=core
#export KMP_AFFINITY=scatter,granularity=core

# could make explicit binding, but machine seems to do the right thing
# export I_MPI_PIN_DOMAIN=numa
mpiexec.hydra -prepend-rank -genvall -n 4 ./dispatch.x $*
