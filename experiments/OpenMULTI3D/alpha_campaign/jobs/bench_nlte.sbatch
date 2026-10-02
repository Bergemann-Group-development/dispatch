#!/bin/bash -l
# Phase-0 timing: 1D NLTE (Mg, Fe) departures plus HRS spectrum (the opacity points
# run in bench_opac.sbatch).
#SBATCH -J alpha_bnl
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node85
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/bench_nlte_opac-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
TIMINGS=$CAMPAIGN/logs/timings_bench_nlte.txt
lscpu | grep -E 'Model name'
cd "$WORK"
for nt in 1 8; do run "$(variant bench/n1d_mg_sun_hrs.nml t$nt)" $nt; done
run "$(variant bench/n1d_mg_giant_hrs.nml t8)" 8
for nt in 8 32; do run "$(variant bench/n1d_fe_sun_hrs.nml t$nt)" $nt; done

echo done
