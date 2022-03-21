#!/bin/bash -l
# --- Log file and error file names ---
#SBATCH -o ./dispatch.log
#SBATCH -e ./dispatch.err
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=debug
#SBATCH --nodes=8                   # Maximum number of nodes to be allocated
#SBATCH --ntasks=8                  # Number of MPI tasks (i.e. processes)
#SBATCH --ntasks-per-node=1         # Maximum number of tasks on each node
#SBATCH --cpus-per-task=32           # Number of cores per MPI task 
# --- User notifications ---
#SBATCH --time 00:30:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=512M
unset I_MPI_PMI_LIBRARY 
mpirun ./dispatch.x