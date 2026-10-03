#!/bin/bash -l
# Cost and accuracy of cheaper 3D stage-2 SED settings: runs val/*.nml (written by
# make_nml.py sed_validate) one after another on all cores, mapping the <3D> departures
# a namelist reads first when val/map.txt lists them.
#   sbatch -w NODE -c 24 --mem=80G --export=ALL,BIN=$C/bin/8ce90476b jobs/sed_validate.sh
#SBATCH -J sed_val
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/sedval-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export TIMINGS=$CAMPAIGN/logs/timings_sedval_$SLURM_JOB_ID.txt
echo "host $(hostname -s), $SLURM_CPUS_PER_TASK threads, BIN=$BIN"
cd $WORK
while read -r nml; do
  [ -z "$nml" ] && continue
  grep -q 'OpenMULTI3D terminated' $CAMPAIGN/logs/$(echo ${nml%.nml} | tr / _).log 2>/dev/null && continue
  src=$(grep -o "atmos_file='[^']*'" "$nml" | cut -d"'" -f2)
  grep " $src " val/map.txt 2>/dev/null | while read -r a b c d; do
    [ -f "$c/save/pop_${d}_meta.txt" ] && continue
    PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/map_departures.py $a $b $c $d
  done
  run "$nml" $SLURM_CPUS_PER_TASK
  tail -1 $TIMINGS
done < val/tasks.txt
