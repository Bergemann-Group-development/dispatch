#!/bin/bash
#SBATCH --job-name=develjob
#SBATCH --partition=astro_devel                                         # devel: 2h, short: 12h, long:120h
#SBATCH --ntasks=16                                                     # number of MPI process
#SBATCH --ntasks-per-node=2 --cpus-per-task=10                          # one per socket
##SBATCH --ntasks-per-node=2 --cpus-per-task=10 --threads-per-core=2     # one per socket + hyperthreading
##SBATCH --ntasks-per-node=1 --cpus-per-task=20                          # one per node

cd $SLURM_SUBMIT_DIR

# time stamp files
run="data"
time=`date +%y%m%d_%H%M`
log=$run/log_$time
x=./dispatch.x
ln -sf $log log

# redirect stdout and stderr
exec >$log 2>&1

# number of threads and OpenMP affinity
export KMP_AFFINITY="granularity=core,scatter"
export OMP_NUM_THREADS=${1:-$SLURM_CPUS_PER_TASK}
echo "SLURM_JOB_ID    = $SLURM_JOB_ID"
echo "KMP_AFFINITY    = $KMP_AFFINITY"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"

# override #threads for tuned hyperthreading
export OMP_NUM_THREADS=15
srun --label --distribution=block:block ./dispatch.x
