#!/bin/bash -l
# Own <3D> averages and their 1D runs, one model after another (claim-based like
# run_queue.sbatch, so jobs on several nodes can share the list).  Per model of
# work/avg3d/models.txt (make_nml.py avg3d): the atmosphere-only snapshot runs,
# avg3d.py, the b4 <3D> LTE, Mg and Fe runs with the departure build (BIN, default
# bin/current as batch 2), then the stage-2 SED with SED_BIN (f63e2fc95 as batch 2).
#   sbatch -w mpia-slurm-node88 -c 1 --mem=4G jobs/avg3d_queue.sbatch
#SBATCH -J alpha_avg3d
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/avg3d-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
SED_BIN=${SED_BIN:-$CAMPAIGN/bin/f63e2fc95}
export TIMINGS=$CAMPAIGN/logs/timings_avg3d_$SLURM_JOB_ID.txt
CLAIMS=${CLAIMS:-$CAMPAIGN/logs/claims_avg3d}
mkdir -p "$CLAIMS"
NT=$SLURM_CPUS_PER_TASK
done_ok () { grep -q 'OpenMULTI3D terminated' "$CAMPAIGN/logs/$(echo ${1%.nml} | tr / _).log" 2>/dev/null; }
cd "$WORK"
while read -r model feh snaps; do
  [ -z "$model" ] && continue
  mkdir "$CLAIMS/$model" 2>/dev/null || continue
  runs=()
  for nml in $snaps; do
    done_ok "$nml" || run "$nml" "$NT"
    runs+=("$CAMPAIGN/runs/${nml%.nml}")
  done
  PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/avg3d.py \
    "$WORK/avg3d/$model.atmos" "$feh" "${runs[@]}" || { echo "$model: average failed"; continue; }
  for el in lte mg fe; do
    [ -f "b4/av3d/$model/$el.nml" ] || continue
    done_ok "b4/av3d/$model/$el.nml" || run "b4/av3d/$model/$el.nml" "$NT"
  done
  [ -f "b4/av3d/$model/fe.nml" ] && python3 $CAMPAIGN/jobs/make_nml.py avg3d_sed "$model" > /dev/null &&
    { done_ok "b4/sed/av3d_$model.nml" || BIN=$SED_BIN run "b4/sed/av3d_$model.nml" "$NT"; }
  echo "$model done $(date +%H:%M)"
done < "$WORK/avg3d/models.txt"
