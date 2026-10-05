#!/bin/bash -l
# Crop the snapshots of crop/tasks.txt (make_nml.py pilot_fe4) with crop3d.py, then
# build any line grid of lgrid/todo.txt that the listed 3D Fe namelists need.
#   sbatch --nice=100 jobs/crop3d.sbatch p1/tasks_3d_fe4.txt
#SBATCH -J alpha_crop
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/crop-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
cd "$WORK"
set -e
while read -r snap run out ltau; do
  [ -f "$out/$(basename $snap)" ] && continue
  PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/crop3d.py "$snap" "$run" "$out" "$ltau"
done < crop/tasks.txt
for nml in $(cat "$1"); do
  g=$(grep -o "line_grid_file='[^']*'" "$nml" | cut -d"'" -f2 | xargs dirname)
  [ -f "$g/line_grid.txt" ] && continue
  cmd=$(grep "^$g " lgrid/todo.txt | cut -d' ' -f2-)
  ( cd $CL && PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/$cmd ) > "$CAMPAIGN/logs/lgrid_$(basename $g).log" 2>&1
done
echo done
