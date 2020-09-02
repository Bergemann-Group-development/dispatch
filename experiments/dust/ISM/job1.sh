#!/bin/bash
#SBATCH --job-name=1pc-2
#SBATCH --exclusive --no-requeue
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=96 --ntasks-per-core=2
##SBATCH --partition=astro2_devel
##SBATCH --partition=astro2_long
#SBATCH --partition=astro2_short
##SBATCH --exclude=astro06-08
##SBATCH --mem=NNNgb

# number of threads and OpenMP affinity
export OMP_NUM_THREADS=48

set -x

# data/run directory
run=tmp
#run=input-copyee2a
#run=1pc-2
#run=1pc-2-128
mkdir -p data/$run

# file names for input and executable
in=$run.nml
exe=./dispatch.x

# time stamp files
time=`date +%y%m%d_%H%M%S`
log=data/$run/log_$time
cp $in data/$run/input_$time.nml

# symbolic link for convenience
ln -sf $log log

# redirect stdout and stderr
exec >$log 2>&1
echo jobid=$SLURM_JOB_ID

#module purge
#module load astro
#module load intel
#module load intelmpi
module list

# No MPI execution command
export KMP_AFFINITY="verbose,granularity=core,scatter"
#srun --cpu_bind=threads $exe $in

# With MPI execution command
srun --cpu_bind=threads --distribution=block:block $exe $in
#srun --cpu_bind=threads --label --distribution=block:block $exe $in
