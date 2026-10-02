#!/bin/bash -l
# Phase-0 timing: LTE opacities at individual (T, rho) points through tumult.hop.
#SBATCH -J alpha_bop
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/bench_opac-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export HOP_LIB=$BIN/libhop.so
ulimit -s unlimited     # ifort keeps large automatic arrays on the main-thread stack
export PYTHONPATH=$CLONE   # tumult from the clone checkout
lscpu | grep -E 'Model name'
OMP_NUM_THREADS=1 /usr/bin/time -f "opacity_bench %e s %M KB rc=%x" \
  "$CLONE/.venv/bin/python" "$CAMPAIGN/jobs/opacity_bench.py"
