#!/bin/bash
#SBATCH --partition=astro_long --ntasks=1 --ntasks-per-node=1
#SBATCH --ntasks-per-core=1 --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --job-name=m05r27p6a

# Don't edit this template file; instead copy it to e.g. long.sh and submit that
# ==============================================================================

# Run name: this & the job-name are the only things that need to be changed, but
# note that the input &io_params namelist must contain a line "datadir='./data'
run=m05r27t10p6a

# -------------------------------- do not edit --------------------------------------------
# The lines below are used to create a track of job input files, executables, and log files
cd $SLURM_SUBMIT_DIR
dir=data/$run
mkdir -p $dir
time=`date +%y%m%d_%H%M`
log=log_$time
in=$run.nml
cp -p $in               $dir/in_$time
cp -p dispatch.x        $dir/dispatch_$time
sed 's:./data:.:' $in > $dir/input.nml

# redirect stdout and stderr
cd $dir
exec >$log 2>&1
echo jobid=$SLURM_JOB_ID
ls -l dispatch_$time

# number of threads and OpenMP affinity
export OMP_NUM_THREADS=${1:-$SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="granularity=core,scatter"
echo $KMP_AFFINITY 
echo "threads=$OMP_NUM_THREADS"
srun --distribution=block:block ./dispatch_$time
# -------------------------------- do not edit --------------------------------------------
