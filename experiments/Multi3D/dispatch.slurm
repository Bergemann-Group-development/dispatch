#!/bin/bash -l
# --- Log file and error file names ---
#SBATCH -o ./dispatch.log
#SBATCH -e ./dispatch.err
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=interactive
#SBATCH --ntasks=1                  # Number of MPI tasks (i.e. processes)
#SBATCH --ntasks-per-node=1         # Maximum number of tasks on each node
#SBATCH --cpus-per-task=8         # Number of cores per MPI task 
#SBATCH --ntasks-per-core=1
# --- User notifications ---
#SBATCH --time 03:00:00
rm ./out/*
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=1024M
unset I_MPI_PMI_LIBRARY 
mpirun ./dispatch.x