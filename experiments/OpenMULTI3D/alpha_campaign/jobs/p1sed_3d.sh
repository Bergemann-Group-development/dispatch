#!/bin/bash -l
# Pilot stage 2 in 3D: refresh the task lists (make_nml.py pilot_sed: every finished 3D
# Mg snapshot, and Mg+Fe where its 3D Fe run has finished), map the <3D> departures onto
# each snapshot (p1sed/map.txt), then the 3D line-list SEDs one after another on all
# cores.  Maps and SEDs that already exist are skipped, so the job can be resubmitted
# as more pilot snapshots finish.  Needs BIN of a build with 7425e2438.
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export TIMINGS=$CAMPAIGN/logs/timings_p1sed3d_$SLURM_JOB_ID.txt
cd $WORK && python3 $CAMPAIGN/jobs/make_nml.py pilot_sed
while read -r a b c d; do
  [ -z "$a" ] && continue
  [ -f "$c/save/pop_${d}_meta.txt" ] && continue
  PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/map_departures.py $a $b $c $d
done < $WORK/p1sed/map.txt
while read -r nml; do
  log=$CAMPAIGN/logs/$(echo ${nml%.nml} | tr / _).log
  grep -q 'OpenMULTI3D terminated' "$log" 2>/dev/null && continue
  run "$nml" $SLURM_CPUS_PER_TASK
done < $WORK/p1sed/tasks_3d.txt
cat $TIMINGS
