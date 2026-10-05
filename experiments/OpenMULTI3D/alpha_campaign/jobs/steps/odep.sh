# STEPS=odep for queue.sbatch: O departures of one MARCS model (make_nml.py o_dep4: line
# grid, [Fe/H]-matched absmet), then the stage-2 SED of the O windows (make_nml.py o_sed).
# A line is "MODEL GRID"; a missing line grid is built first from lgrid/todo.txt. BIN and
# SED_BIN default to the builds of the b6 Mg/Fe runs (steps/dep4.sh).
#   sbatch -w NODE -c 4 --mem=8G --export=ALL,STEPS=odep,BIN=$C/bin/23675093a \
#     jobs/queue.sbatch b7/tasks_dep.txt 1
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
SED_BIN=${SED_BIN:-$CAMPAIGN/bin/f63e2fc95}

task_key () { echo "$1"; }
# done: the SED terminated and read lines (a window file the reader mishandles gives
# "0 of N lines read" and a line-free spectrum)
task_done () {
  local log=$CAMPAIGN/logs/b7_sed_$1.log
  grep -qs 'OpenMULTI3D terminated' "$log" && ! grep -qs '^0 of .* lines read from linelist .*vald' "$log"
}
task () {
  local model=$1 g=$2 cmd
  if [ ! -f "$g/line_grid.txt" ]; then
    cmd=$(grep "^$g " lgrid/todo.txt | cut -d' ' -f2-)
    ( cd $CL && PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/$cmd ) \
      > "$CAMPAIGN/logs/lgrid_$(basename $g).log" 2>&1 || { echo "$model: line grid $g failed"; return 1; }
  fi
  grep -qs 'OpenMULTI3D terminated' "$CAMPAIGN/logs/b7_n1d_${model}_o.log" || run "b7/n1d/${model}_o.nml" "$NT"
  python3 $CAMPAIGN/jobs/make_nml.py o_sed "$model" > /dev/null && BIN=$SED_BIN run "b7/sed/$model.nml" "$NT"
  task_done "$model"
}
