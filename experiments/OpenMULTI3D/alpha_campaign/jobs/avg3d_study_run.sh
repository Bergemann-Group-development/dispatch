#!/bin/bash -l
# avg3d_study_run.sh LIST THREADS: run the namelists of LIST (relative to work/) one after
# another with THREADS OpenMP threads, skipping finished ones and ones another runner has
# claimed (atomic mkdir under logs/claims_avg3d_study).  Works under Slurm and on the vis
# nodes outside it, so several machines can share one list:
#   sbatch -w mpia-slurm-node86 --nice -c 32 --mem=120G -t 2-0 \
#     --export=ALL,BIN=$C/bin/4890f9fde jobs/avg3d_study_run.sh avg3d_study/tasks_3d.txt 32
#   ssh vis2 "cd $C && BIN=$C/bin/4890f9fde nohup jobs/avg3d_study_run.sh LIST 36 > logs/avg3d_study_vis2.out 2>&1 &"
# A claim left by a runner that died is removed by hand (rmdir).
#SBATCH -J avg3d_study
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/avg3d_study-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
LIST=$1
NT=$2
CLAIMS=${CLAIMS:-$CAMPAIGN/logs/claims_avg3d_study}
export TIMINGS=$CAMPAIGN/logs/timings_avg3d_study_$(hostname -s).txt
mkdir -p "$CLAIMS"
cp "$WORK/$LIST" "$CLAIMS/.list.$$"
while read -r nml; do
  [ -z "$nml" ] && continue
  key=${nml%.nml}; key=${key//\//_}
  grep -qs 'OpenMULTI3D terminated' "$CAMPAIGN/logs/$key.log" && continue
  mkdir "$CLAIMS/$key" 2>/dev/null || continue
  echo "$(hostname -s) ${SLURM_JOB_ID:-$$}" > "$CLAIMS/$key/owner"
  echo "$(date +%H:%M) start $nml"
  run "$nml" "$NT"
  grep "^${nml%.nml} " "$TIMINGS" | tail -1
done < "$CLAIMS/.list.$$"
rm -f "$CLAIMS/.list.$$"
