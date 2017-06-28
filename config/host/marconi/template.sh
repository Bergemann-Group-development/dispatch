#!/bin/bash
#PBS -l select=4:ncpus=68:mpiprocs=1:mem=64GB
#PBS -m b -M aake@nbi.ku.dk
#PBS -l walltime=0:30:00
#PBS -A Pra13_3286

cd $PBS_JOBDIR

# time stamp log file: data/log_yymmdd_HHMM
run="data"
time=`date +%y%m%d_%H%M`
log=$run/log_$time
x=./dispatch.x
ln -sf $log log

# redirect stdout and stderr
exec >$log 2>&1

# number of threads and OpenMP affinity
export KMP_AFFINITY="granularity=core,scatter"
export OMP_NUM_THREADS=272

# execute with numactl memory placement
mpirun -np 4 numactl ./dispatch.x
