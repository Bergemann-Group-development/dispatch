#!/bin/bash
# Batch 2: the paired pilot on seven Stagger models with K snapshots each (default 2):
# 1D LTE/NLTE on the MARCS twins and <3D> averages, 3D LTE, 3D NLTE Mg and 3D NLTE Fe.
# Run on the login node: it only writes namelists and submits.
#
# The first 3D Mg run (the Sun, first snapshot) runs on its own; a gate job then checks
# it (check_runs.py --gate: termination, convergence, finite positive departures, a sane
# NLTE-LTE difference) and the 3D Fe array starts only if the gate passes (afterok).
# A failed gate leaves the Fe array pending with DependencyNeverSatisfied; resubmit it by
# hand after a look at logs/gate-*.out.  Walltimes allow about twice the expected time.
set -e
K=${1:-2}
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
cd "$CAMPAIGN/work"
python3 ../jobs/make_nml.py pilot "$K"
head -n 1 p1/tasks_3d_mg.txt > p1/tasks_3d_mg_first.txt
tail -n +2 p1/tasks_3d_mg.txt > p1/tasks_3d_mg_rest.txt
first=$(head -n 1 p1/tasks_3d_mg_first.txt)
n3d=$(wc -l < p1/tasks_3d_fe.txt)
cd "$CAMPAIGN"
sub () { sbatch --parsable "$@"; }
# whole nodes, as queue.sbatch defaults to one core
node="--exclusive -c 32 --mem=110G"
# 42 1D runs, 8 threads each, 4 at a time: ~10 min
j1d=$(sub -t 01:00:00 -J p1_1d $node jobs/queue.sbatch p1/tasks_1d.txt 8)
# the gate's 3D Mg run: ~2-3 h (convlim 1e-4, up to 150 iterations of ~75 s, then spectra)
jmg0=$(sub -t 08:00:00 -J p1_3dmg0 $node jobs/queue.sbatch p1/tasks_3d_mg_first.txt 32)
# 3D LTE, one run after the other at 32 threads, ~15 min each at R=2e5
jlte=$(sub -t 10:00:00 -J p1_3dlte $node jobs/queue.sbatch p1/tasks_3d_lte.txt 32)
# the other 3D Mg runs in two shards
jmg=$(sub --array=0-1 -t 40:00:00 -J p1_3dmg $node jobs/queue.sbatch p1/tasks_3d_mg_rest.txt 32)
jgate=$(sub --dependency=afterany:$jmg0 -J p1_gate --partition=four-wks --nodes=1 \
  --nodelist='mpia-slurm-node[81-90]' --cpus-per-task=2 --mem=16G -t 00:30:00 \
  --output="$CAMPAIGN/logs/gate-%j.out" \
  --wrap="PYTHONPATH=$CLONE OMP_NUM_THREADS=1 $CLONE/.venv/bin/python $CAMPAIGN/jobs/check_runs.py --gate $CAMPAIGN/runs/${first%.nml}")
# 3D NLTE Fe, one run per array element (~80 GB, one per node), first snapshots first
jfe=$(sub --dependency=afterok:$jgate --array=0-$((n3d - 1)) -t 48:00:00 -J p1_3dfe \
  $node jobs/queue.sbatch p1/tasks_3d_fe.txt 32)
echo "1d $j1d, 3dmg0 $jmg0, 3dlte $jlte, 3dmg $jmg, gate $jgate, 3dfe $jfe" | tee -a logs/submitted_pilot.txt
