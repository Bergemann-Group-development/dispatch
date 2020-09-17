#!/bin/bash

# Pick up run name from argument, and set time string
run=$1
time=`date +%y%m%d_%H%M%S`

# Set queue to short by default, or 2nd argument
q=${2:-short}

# File names for input, log, and executable
in=$run.nml
log=data/$run/log_$time
exe=./dispatch.x

# Make data/run directory
mkdir -p data/$run
cp $in data/$run/input_$time.nml
ln -sf $log log

# Submit job script
sbatch << EOF
#!/bin/bash
#SBATCH --job-name=$run
#SBATCH --exclusive --no-requeue
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=96 --ntasks-per-core=2
#SBATCH --partition=astro2_$q
##SBATCH --exclude=astro06-08
##SBATCH --mem=NNNgb

set -x

# number of threads and OpenMP affinity
export OMP_NUM_THREADS=48
export KMP_AFFINITY="verbose,granularity=core,scatter"

# redirect stdout and stderr
exec >$log 2>&1
echo jobid=\$SLURM_JOB_ID

#module purge
#module load astro
#module load intel
#module load intelmpi
module list

# No MPI execution command
srun --cpu_bind=threads $exe $in

# With MPI execution command
#srun --cpu_bind=threads --distribution=block:block $exe $in
#srun --cpu_bind=threads --label --distribution=block:block $exe $in
EOF
