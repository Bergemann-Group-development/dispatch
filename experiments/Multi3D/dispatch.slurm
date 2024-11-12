#!/bin/bash -l
# --- Log file and error file names ---
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH --output ./slurm/slurm-%j.out
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=four-wks
#SBATCH --nodes=1                  # Number of MPI tasks (i.e. processes)
##SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --exclude=mpia-slurm-node[71,74]
#SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node
#SBATCH --cpus-per-task=32         # Number of cores per MPI task 
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive
#SBATCH --mem=100G
# --- User notifications ---
##SBATCH --time 24:00:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_STACKSIZE=256M
#unset OMP_NUM_THREADS
unset I_MPI_PMI_LIBRARY

for arg in "$@"
do
  mpirun ./dispatch.x $arg
done