#!/bin/bash -l
# Queue runner: take the runs of a task list (namelists relative to work/) one after
# another, each on all the job's cores, until the list is exhausted.  Several jobs can
# serve one list: a run is claimed by an atomic mkdir under logs/claims/, and finished or
# claimed runs are skipped.  Size the job with sbatch options, e.g. the four free cores
# of node83:
#   sbatch -w mpia-slurm-node83 -c 4 --mem=30G --nice=100 jobs/run_queue.sbatch x3d/tasks.txt
# Node83's /opt/intel is a stale mount, so the job puts the campaign's copy of the
# OpenMP runtime (lib/libiomp5.so, from the login node's oneAPI 2024.0) on the path.
#SBATCH -J alpha_queue
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/queue-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
LIST=$1
export TIMINGS=$CAMPAIGN/logs/timings_$(basename "$(dirname "$WORK/$LIST")")_queue_$SLURM_JOB_ID.txt
mkdir -p "$CAMPAIGN/logs/claims"
echo "host $(hostname -s), $SLURM_CPUS_PER_TASK threads per run"
while read -r nml; do
  [ -z "$nml" ] && continue
  name=${nml%.nml}
  mkdir "$CAMPAIGN/logs/claims/${name//\//_}" 2>/dev/null || continue
  run "$nml" "$SLURM_CPUS_PER_TASK"
  tail -1 "$TIMINGS"
done < "$WORK/$LIST"
