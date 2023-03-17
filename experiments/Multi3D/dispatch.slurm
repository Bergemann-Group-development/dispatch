#!/bin/bash -l
# --- Log file and error file names ---
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=interactive
#SBATCH --nodes=1                  # Number of MPI tasks (i.e. processes)
#SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node
#SBATCH --cpus-per-task=30         # Number of cores per MPI task 
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive
# --- User notifications ---
#SBATCH --time 3-00:00:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_STACKSIZE=256M
#unset OMP_NUM_THREADS
unset I_MPI_PMI_LIBRARY
mpirun ./dispatch.x $1
#mpirun ./dispatch.x input_fe_nlte.nml