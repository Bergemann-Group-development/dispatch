#!/bin/bash -l
# delta_maps.py (or pilot_compare.py for mode 'pilot') for the given modes into runs/delta/<MODE><SUFFIX>.csv, e.g.
#   sbatch --nice=100 jobs/delta_maps.sbatch A2 N W X
#   sbatch --nice=100 --export=ALL,LTE_IT1=1,SUFFIX=_it1 jobs/delta_maps.sbatch N W
# The environment (AV3D, AVG_DIR, B4_NLTE, LTE_IT1, OLD_SED, ONLY) passes through to the scripts.
#SBATCH -J alpha_delta
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/delta-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
mkdir -p $CAMPAIGN/runs/delta
cd $CAMPAIGN/jobs
for m in "$@"; do
  script=delta_maps.py; arg=$m
  [ "$m" = pilot ] && { script=pilot_compare.py; arg=; }     # mode 'pilot': pilot_compare.py
  PYTHONPATH=$CL:$CAMPAIGN/jobs OMP_NUM_THREADS=1 $CL/.venv/bin/python $script $arg \
    $CAMPAIGN/runs/delta/${m}${SUFFIX:-_p5}.csv
done
