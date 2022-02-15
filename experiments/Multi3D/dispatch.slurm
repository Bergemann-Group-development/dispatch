#!/bin/bash -l
# --- Log file and error file names ---
#SBATCH -o ./dispatch.log
#SBATCH -e ./dispatch.err
# --- Job name and working directory ---
#SBATCH -D ./
#SBATCH -J M3D_DIS
# --- CPU and node logistics ---
#SBATCH --partition=astro_devel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem=8000
# --- User notifications ---
#SBATCH --time 00:05:00
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -n 1 dispatch.x