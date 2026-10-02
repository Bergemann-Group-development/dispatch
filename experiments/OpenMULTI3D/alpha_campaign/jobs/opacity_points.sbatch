#!/bin/bash -l
# LTE opacity points (batch 1b): sample the (T, rho) domain once, then compute every
# (composition, shard) pair as independent one-thread processes filling the node.
#   sbatch jobs/opacity_points.sbatch "-4 -3 -2 -1.5 -1 -0.5 0 0.5" 64
#SBATCH -J alpha_opac
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=110G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/opac-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export PYTHONPATH=$CLONE HOP_LIB=$BIN/libhop.so OMP_NUM_THREADS=1
ulimit -s unlimited
PY="$CLONE/.venv/bin/python $CAMPAIGN/jobs/opacity_points.py"
FEHS=${1:-"-4 -3 -2 -1.5 -1 -0.5 0 0.5"}; NSHARD=${2:-64}   # existing shards are skipped
[ -f "$CAMPAIGN/runs/opac/points.npz" ] || $PY points || exit 1
for f in $FEHS; do for k in $(seq 0 $((NSHARD-1))); do echo "$f $k $NSHARD"; done; done | \
  xargs -P "$SLURM_CPUS_PER_TASK" -L 1 /usr/bin/time -f "%e s %M KB rc=%x" $PY shard
ls "$CAMPAIGN"/runs/opac/feh*/ | grep -c npz; du -sh "$CAMPAIGN/runs/opac"
