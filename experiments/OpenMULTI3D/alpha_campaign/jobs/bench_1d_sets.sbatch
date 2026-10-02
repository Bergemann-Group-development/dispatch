#!/bin/bash -l
# Phase-0 timing: marginal cost of abundance-offset sets (set_file) in one 1D LTE run,
# and the run without linelist_out.txt.
#SBATCH -J alpha_bset
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node82
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/bench_1d_sets-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
TIMINGS=$CAMPAIGN/logs/timings_bench_1d.txt
cd "$WORK"
run "$(variant bench/l1d_sun_hrs_sets.nml t1)" 1
run "$(variant bench/l1d_sun_hrs_nolinelist_out.nml t1)" 1
run "$(variant bench/l1d_sun_hrs.nml t1rep)" 1
