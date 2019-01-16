#!/bin/bash

# Run name: this & the job-name are the only things that need to be changed, but
# note that the input &io_params namelist must contain a line "datadir='./data'
run=$1
#run=1024m16_hllc

# -------------------------------- do not edit --------------------------------------------
# The lines below are used to create a track of job input files, executables, and log files
cd $SLURM_SUBMIT_DIR
dir=data/$run
mkdir -p $dir
ln -sf $dir run
time=`date +%y%m%d_%H%M`
log=$dir/log_$time
in=$run.nml
x=${2:-dispatch.x}
#x=dispatch_hllc.x
cp -p $in               $dir/in_$time
cp -p $x                $dir/dispatch_$time
ln -sf log_$time $dir/log
mv prv pprv || echo "no pprv"
mv log prv || echo "no prv"
ln -sf $log log

# start a jobhelper in the data dir, unless a process is already running
#ps -u $USER | grep jobhelper || ( cd $dir; jobhelper >jobhelper.log 2>&1 </dev/null & )
jobhelper_start $dir

# redirect stdout and stderr
exec >$log 2>&1
echo jobid=$SLURM_JOB_ID
ls -l $x $dir/dispatch_$time

# number of threads and OpenMP affinity
export OMP_NUM_THREADS=${3:-$SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="granularity=core,scatter"
export KMP_STACKSIZE="30m"
echo "KMP_AFFINITY=$KMP_AFFINITY"
echo "threads=$OMP_NUM_THREADS"
#srun --label --distribution=block:block $dir/dispatch_$time $in
srun --distribution=block:block $dir/dispatch_$time $in

# ad hoc
#export OMP_NUM_THREADS=20
#srun -n 64 -c 20 --distribution=block:block $dir/dispatch_$time $in
mv $dir ${dir}_$time
# -------------------------------- do not edit --------------------------------------------
