#!/bin/bash -l
# Phase-0 timing: 3D LTE HRS synthesis and 3D NLTE Mg departures on one Stagger-grid
# snapshot (t5777g44m00, 80x80x240 interpolated to the 3d preset's 30x30x96).
#SBATCH -J alpha_b3d
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node86
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/bench_3d-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
TIMINGS=$CAMPAIGN/logs/timings_bench_3d.txt
lscpu | grep -E 'Model name'
cd "$WORK"
run "$(variant bench/l3d_sun_hrs_R1e5.nml t32)" 32
run "$(variant bench/n3d_mg_sun_dep.nml t32)" 32
echo done
