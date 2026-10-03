#!/bin/bash -l
# The vmic-grid 1D LTE partners of the 3D LTE snapshots (make_nml.py vmic_grid, run
# here first so that it sees the 3D LTE runs and own <3D> averages finished by then),
# claimed one by one like run_queue.sbatch.  Submit after the avg3d queue jobs.
#SBATCH -J alpha_vmic
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/vmic-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TIMINGS=$CAMPAIGN/logs/timings_vmic_$SLURM_JOB_ID.txt
CLAIMS=$CAMPAIGN/logs/claims_vmic
mkdir -p "$CLAIMS"
cd "$WORK"
[ -n "$REGEN" ] && python3 $CAMPAIGN/jobs/make_nml.py vmic_grid
while read -r nml; do
  [ -z "$nml" ] && continue
  grep -q 'OpenMULTI3D terminated' "$CAMPAIGN/logs/$(echo ${nml%.nml} | tr / _).log" 2>/dev/null && continue
  mkdir "$CLAIMS/$(echo ${nml%.nml} | tr / _)" 2>/dev/null || continue
  run "$nml" "$SLURM_CPUS_PER_TASK"
done < "$WORK/b5/tasks_vmic.txt"
