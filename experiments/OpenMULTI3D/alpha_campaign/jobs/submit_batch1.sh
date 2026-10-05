#!/bin/bash
# Batch 1 (Friday night): 1a, the 1D LTE grid over the 4MOST domain, and 1b, LTE opacity
# points over the (T, rho) domain of the atmospheres.  Run on the login node: it only
# writes namelists and submits.  Expected: 1a ~15 min on three nodes (2354 runs,
# ~21 core-h), 1b ~15 min on one node (65536 points, ~50 GB).
set -e
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
cd "$CAMPAIGN/work"
python3 ../jobs/make_nml.py batch1
cd "$CAMPAIGN"
sbatch --array=0-2 --exclusive -c 32 --mem=110G -t 01:00:00 -J b1_1dlte jobs/queue.sbatch b1/tasks.txt 1
sbatch -t 01:00:00 -J b1_opac jobs/opacity_points.sbatch "-4 -3 -2 -1.5 -1 -0.5 0 0.5" 64
