#!/bin/bash -l
# Sanity checks of finished runs (check_runs.py), 16 runs per process, one process per
# core.  RUNS is a task list of namelists relative to work/ (each run's folder is
# runs/<path without .nml>).  Output: logs/check-<job>.out, one line per run.
#   sbatch jobs/check_runs.sbatch b1/tasks.txt
#SBATCH -J alpha_check
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/check-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export PYTHONPATH=$CLONE OMP_NUM_THREADS=1
sed "s#^#$CAMPAIGN/runs/#; s#\.nml\$##" "$WORK/$1" | \
  xargs -P "$SLURM_CPUS_PER_TASK" -n 16 "$CLONE/.venv/bin/python" "$CAMPAIGN/jobs/check_runs.py" \
  | sort > "$CAMPAIGN/logs/check-$SLURM_JOB_ID.txt"
echo "$(grep -c '^ok' "$CAMPAIGN/logs/check-$SLURM_JOB_ID.txt") ok, $(grep -c '^FAIL' "$CAMPAIGN/logs/check-$SLURM_JOB_ID.txt") failed"
grep '^FAIL' "$CAMPAIGN/logs/check-$SLURM_JOB_ID.txt" | head -50
