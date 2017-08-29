#!/bin/bash

cd $SLURM_SUBMIT_DIR
set -x
env > env.log

# file names for input and executable
in=$1
bin=./
exe=dispatch.x
touch restart.flag

# time stamp files
run=data/`basename $in .nml`
time=`date +%y%m%d_%H%M`
log=$run/log_$time
x=$run/${exe}_$time
mv log prv || echo "no prv"
ln -sf $log log
cp -p $in $run/in_$time
cp -p $bin/$exe $x

# redirect stdout and stderr
exec >$log 2>&1
echo jobid=$SLURM_JOB_ID
ls -lL $bin/$exe $log

# number of threads and OpenMP affinity
export OMP_NUM_THREADS=${2:-$SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="granularity=core,scatter,1,0"
echo OMP_NUM_THREADS=$OMP_NUM_THREADS

#export I_MPI_DEBUG="3"
#srun --label --distribution=block:block $x $in
srun --distribution=block:block $x $in
