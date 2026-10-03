#!/bin/bash -l
# Pilot stage 2 in 3D: refresh the task lists (make_nml.py pilot_sed [RT_MAXITER]: every
# finished 3D Mg snapshot, and Mg+Fe where its 3D Fe run has finished, in priority
# order), then the 3D line-list SEDs one after another on all cores, each preceded by
# the <3D> departure maps it reads (p1sed/map.txt).  A SED is claimed by an atomic mkdir
# under logs/claims_p1sed/, so several jobs on different nodes can share the list;
# finished SEDs are skipped.  Needs BIN of a build with 7425e2438 (8ce90476b reads one
# Jnu set per final ray).  With LOCAL=1 a SED runs on the node's /scratch when it has
# 150 GB free (its Jnu_out files reach about 120 GB; beegfs throughput varies with the
# other users' load); everything but the Jnu files is copied back to runs/ afterwards.
#   sbatch -w NODE -c 24 --mem=80G --time=4-00:00:00 --export=ALL,BIN=$C/bin/8ce90476b jobs/p1sed_3d.sh
#SBATCH -J p1sed_3d
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/p1sed3d-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export TIMINGS=$CAMPAIGN/logs/timings_p1sed3d_$SLURM_JOB_ID.txt
CLAIMS=${CLAIMS:-$CAMPAIGN/logs/claims_p1sed}
mkdir -p "$CLAIMS"
echo "host $(hostname -s), $SLURM_CPUS_PER_TASK threads, BIN=$BIN"
# the output scattering iterations: RT_MAXITER, else work/p1sed/rt_maxiter, else 3
[ -z "$RT_MAXITER" ] && [ -f $WORK/p1sed/rt_maxiter ] && RT_MAXITER=$(cat $WORK/p1sed/rt_maxiter)
cd $WORK && python3 $CAMPAIGN/jobs/make_nml.py pilot_sed ${RT_MAXITER:-3} || { echo 'pilot_sed failed; stale task list, stopping'; exit 1; }
while read -r nml; do
  [ -z "$nml" ] && continue
  log=$CAMPAIGN/logs/$(echo ${nml%.nml} | tr / _).log
  grep -q 'OpenMULTI3D terminated' "$log" 2>/dev/null && continue
  mkdir "$CLAIMS/$(echo ${nml%.nml} | tr / _)" 2>/dev/null || continue
  src=$(grep -o "atmos_file='[^']*'" "$nml" | cut -d"'" -f2)
  grep " $src " $WORK/p1sed/map.txt | while read -r a b c d; do
    [ -f "$c/save/pop_${d}_meta.txt" ] && continue
    PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/map_departures.py $a $b $c $d
  done
  free=$(df -BG --output=avail /scratch 2>/dev/null | tail -1 | tr -dc 0-9)
  if [ "${LOCAL:-0}" = 1 ] && [ "${free:-0}" -ge 150 ]; then
    tmp=/scratch/hoppe_alpha_$SLURM_JOB_ID
    mkdir -p $tmp/work/$(dirname $nml)
    ln -sfn $WORK/input $tmp/work/input; ln -sfn $WORK/windows $tmp/work/windows
    sed "s#datadir='[^']*'#datadir='$tmp/runs'#" $nml > $tmp/work/$nml
    ( WORK=$tmp/work run "$nml" $SLURM_CPUS_PER_TASK )
    mkdir -p $CAMPAIGN/runs/${nml%.nml}
    rsync -a --exclude='Jnu*' --exclude='exit_plane*' $tmp/runs/${nml%.nml}/ $CAMPAIGN/runs/${nml%.nml}/
    rm -rf $tmp
  else
    run "$nml" $SLURM_CPUS_PER_TASK
  fi
  tail -1 $TIMINGS
done < $WORK/p1sed/tasks_3d.txt
