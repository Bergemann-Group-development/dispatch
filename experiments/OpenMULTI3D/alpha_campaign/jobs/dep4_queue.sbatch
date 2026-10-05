#!/bin/bash -l
# Mg and Fe departures of the own <3D> averages with the phase-4 setup (make_nml.py
# dep4: line-grid quadrature, [Fe/H]-matched absmet), then the stage-2 SED, one model
# after another; models are claimed like run_queue.sbatch.  Missing line grids are
# built first (lgrid/todo.txt).  BIN needs the big-endian absmet fix (d2f3eabfb), SED_BIN
# defaults to the batch-2 stage-2 build.  Runs of the first setup were moved to
# runs/b4_defq/ before this queue started.
#   sbatch -w NODE -c 1 --mem=4G --export=ALL,BIN=$C/bin/23675093a jobs/dep4_queue.sbatch
# KIND=marcs runs the batch-2 MARCS redo instead (make_nml.py batch2_dep4: b6/n1d and
# b6/sed, list b6/tasks_dep.txt); LIST names another model list of the same kind
# (make_nml.py dep4_list).  An array of 1-core jobs shares a list through the claims.
#SBATCH -J alpha_dep4
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/dep4-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
SED_BIN=${SED_BIN:-$CAMPAIGN/bin/f63e2fc95}
export TIMINGS=$CAMPAIGN/logs/timings_dep4_${SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID}.txt
KIND=${KIND:-av3d}
if [ "$KIND" = marcs ]; then LIST=${LIST:-b6/tasks_dep.txt}; else LIST=${LIST:-b4/tasks_dep.txt}; fi
CLAIMS=$CAMPAIGN/logs/claims_dep4$([ "$KIND" = marcs ] && echo _b6)
mkdir -p "$CLAIMS"
NT=$SLURM_CPUS_PER_TASK
echo "host $(hostname -s), BIN=$BIN"
cd "$WORK"
while read -r model grids; do
  [ -z "$model" ] && continue
  mkdir "$CLAIMS/$model" 2>/dev/null || continue
  ok=1
  for g in $grids; do
    [ -f "$g/line_grid.txt" ] && continue
    cmd=$(grep "^$g " lgrid/todo.txt | cut -d' ' -f2-)
    ( cd $CL && PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/$cmd ) > "$CAMPAIGN/logs/lgrid_$(basename $g).log" 2>&1 \
      || { echo "$model: line grid $g failed"; ok=0; }
  done
  [ $ok = 1 ] || continue
  if [ "$KIND" = marcs ]; then
    for el in mg fe; do run "b6/n1d/${model}_$el.nml" "$NT"; done
    python3 $CAMPAIGN/jobs/make_nml.py b6_sed "$model" > /dev/null && BIN=$SED_BIN run "b6/sed/$model.nml" "$NT"
  else
    for el in mg fe; do run "b4/av3d/$model/$el.nml" "$NT"; done
    python3 $CAMPAIGN/jobs/make_nml.py avg3d_sed "$model" > /dev/null && BIN=$SED_BIN run "b4/sed/av3d_$model.nml" "$NT"
  fi
  echo "$model done $(date +%H:%M)"
done < "$WORK/$LIST"
