#!/bin/bash -l
# Phase-0 timing: 1D LTE and <3D> LTE syntheses of the 4MOST HRS windows, thread
# scaling, line-list variants, hop_prune.x, and a packed node (32 one-thread runs).
#SBATCH -J alpha_b1d
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node82
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=01:30:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/bench_1d-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
TIMINGS=$CAMPAIGN/logs/timings_bench_1d.txt
lscpu | grep -E 'Model name'
cd "$WORK"

# line-list pruning for the sun, then the pruned-list run
OMP_NUM_THREADS=1 /usr/bin/time -f "prune_sun_hrs 1 %e %M %x $(hostname -s)" -a -o "$TIMINGS" \
  "$BIN/hop_prune.x" bench/prune_sun_hrs.nml > "$CAMPAIGN/logs/bench_prune_sun_hrs.log" 2>&1

for nt in 1 4 8 32; do run "$(variant bench/l1d_sun_hrs.nml t$nt)" $nt; done
for n in l1d_sun_hrs_full l1d_sun_hrs_R1e5 l1d_giant_hrs l1d_mpdwarf_hrs l1d_sun_lrs \
         l1d_sun_hrs_pruned av3d_sun_hrs; do
  run "$(variant bench/$n.nml t1)" 1
done
run "$(variant bench/l1d_giant_hrs.nml t8)" 8

# packed node: 32 concurrent one-thread runs (alternating sun and giant)
for i in $(seq -w 1 32); do
  src=bench/l1d_sun_hrs.nml; [ $((10#$i % 2)) -eq 0 ] && src=bench/l1d_giant_hrs.nml
  nml=$(variant $src pack$i)
  run "$nml" 1 &
done
wait
echo done
