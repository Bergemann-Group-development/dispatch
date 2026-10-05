# STEPS=dep4 for queue.sbatch: Mg and Fe departures of one model with the phase-4 setup
# (make_nml.py dep4: line-grid quadrature, [Fe/H]-matched absmet), then its stage-2 SED.
# A line is "MODEL GRID..."; missing line grids are built first from lgrid/todo.txt.
# BIN needs the big-endian absmet fix (d2f3eabfb); SED_BIN defaults to the batch-2
# stage-2 build. KIND=marcs runs the batch-2 MARCS redo (make_nml.py batch2_dep4: b6/n1d
# and b6/sed, list b6/tasks_dep.txt) instead of the own <3D> averages (b4, list
# b4/tasks_dep.txt, make_nml.py dep4_list).
#   sbatch -w NODE -c 1 --mem=4G --export=ALL,STEPS=dep4,BIN=$C/bin/23675093a \
#     jobs/queue.sbatch b4/tasks_dep.txt
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
SED_BIN=${SED_BIN:-$CAMPAIGN/bin/f63e2fc95}
KIND=${KIND:-av3d}

task_key () { echo "$1"; }
task_done () {
  if [ "$KIND" = marcs ]; then local sed=b6/sed/$1; else local sed=b4/sed/av3d_$1; fi
  grep -qs 'OpenMULTI3D terminated' "$CAMPAIGN/logs/${sed//\//_}.log"
}
task () {
  local model=$1 g cmd el; shift
  for g in "$@"; do
    [ -f "$g/line_grid.txt" ] && continue
    cmd=$(grep "^$g " lgrid/todo.txt | cut -d' ' -f2-)
    ( cd $CL && PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/$cmd ) \
      > "$CAMPAIGN/logs/lgrid_$(basename $g).log" 2>&1 || { echo "$model: line grid $g failed"; return 1; }
  done
  if [ "$KIND" = marcs ]; then
    for el in mg fe; do run "b6/n1d/${model}_$el.nml" "$NT"; done
    python3 $CAMPAIGN/jobs/make_nml.py b6_sed "$model" > /dev/null && BIN=$SED_BIN run "b6/sed/$model.nml" "$NT"
  else
    for el in mg fe; do run "b4/av3d/$model/$el.nml" "$NT"; done
    python3 $CAMPAIGN/jobs/make_nml.py avg3d_sed "$model" > /dev/null && BIN=$SED_BIN run "b4/sed/av3d_$model.nml" "$NT"
  fi
  echo "$model done $(date +%H:%M)"
  task_done "$model"
}
