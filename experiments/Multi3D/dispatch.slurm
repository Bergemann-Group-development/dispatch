#!/bin/bash -l
# --- Log file and error file names ---
#SBATCH -o ./dispatch.log
#SBATCH -e ./dispatch.err
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=astro2_devel
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# --- User notifications ---
#SBATCH --time 00:20:00
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=512m
srun dispatch.x