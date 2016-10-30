#!/bin/bash
#PBS -N rt_speed
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:02:00             

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

echo 'job id ' $PBD_JOBID
export OMP_NUM_THREADS=1

# aprun arguments:
# -n: total no of PEs (Processing Elements; MPI tasks)
# -N: number of PEs per node
# -d: number of OpenMP threads per node
# -j: no. PEs per core (hyper-threading)
# -cc: CPU affinity control (choices: cpu, numa_node, none, custom [0,4,3,2,etc.])
# -r:  core specialisation; reserve a number of cores, x, for system/service work
# -ss: memory affinity control; a PE can only allocate memory to a specific "numa_node"

# Launch the parallel job to the allocated compute nodes
make clean
make COMPILER=ifort HW="-qopt-report -vec-threshold0" test.x
aprun -n $PBS_NUM_NODES -N 1 -d $OMP_NUM_THREADS -j 1 ./test.x | tee -a performance/hazelhen.log
