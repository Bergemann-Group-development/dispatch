#!/bin/bash -l
# Array runner: shard k of N takes every N-th line of a task list (namelists relative to
# work/) and runs them PER_TASK_THREADS-threaded, filling the node's cores.
#   sbatch --array=0-3 jobs/run_tasks.sbatch b1/tasks.txt 1
# The task list is sorted so each shard gets the same mix of short and long runs.
#SBATCH -J alpha_tasks
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=110G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/tasks-%A_%a.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
LIST=$1; NT=${2:-1}
N=${SLURM_ARRAY_TASK_COUNT:-1}; K=${SLURM_ARRAY_TASK_ID:-0}
export TIMINGS=$CAMPAIGN/logs/timings_$(basename "$(dirname "$WORK/$LIST")")_${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}_$K.txt
P=$(( SLURM_CPUS_PER_TASK / NT ))
echo "host $(hostname -s) shard $K/$N, $P concurrent runs x $NT threads"
cd "$WORK"
awk -v n="$N" -v k="$K" '(NR-1) % n == k' "$LIST" | \
  xargs -P "$P" -I{} "$CAMPAIGN/jobs/run_one.sh" {} "$NT"
echo "done: $(wc -l < "$TIMINGS") runs, $(awk '$5!=0' "$TIMINGS" | wc -l) nonzero exits"
