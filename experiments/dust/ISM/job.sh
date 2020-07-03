#!/bin/bash
#SBATCH --mem-per-cpu=3g --exclusive
#SBATCH --job-name=dust_ism
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=20
##SBATCH --ntasks-per-core=2
#SBATCH --partition=astro_long
##SBATCH --time=4-23:59:59
##SBATCH --partition=astro_short
#SBATCH --time=10:00:00
#SBATCH --no-requeue
##SBATCH --exclude=astro06-08
##SBATCH --mem=NNNgb

# number of threads and OpenMP affinity
export KMP_AFFINITY="granularity=thread,scatter,1,0"

set -x

# file names for input and executable
exe=./dispatch.x

in=1pc-42.nml

# time stamp files
time=`date +%y%m%d_%H%M`
log=log_$time

# redirect stdout and stderr
exec >$log 2>&1
echo jobid=$SLURM_JOB_ID

#module purge
#module load astro
#module load intel
#module load intelmpi
module list

# No MPI execution command
export OMP_NUM_THREADS=20
srun --cpu_bind=threads $exe $in

# With MPI execution command
#srun --cpu_bind=threads --label --distribution=block:block $exe $in
