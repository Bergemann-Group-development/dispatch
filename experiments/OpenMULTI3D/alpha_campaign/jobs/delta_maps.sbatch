#!/bin/bash -l
# delta_maps.py for the given modes into runs/delta/<MODE><SUFFIX>.csv, e.g.
#   sbatch --nice=100 jobs/delta_maps.sbatch A2 N W X
#   sbatch --nice=100 --export=ALL,LTE_IT1=1,SUFFIX=_it1 jobs/delta_maps.sbatch N W
# The environment (B4_NLTE, LTE_IT1, OLD_SED, ONLY) passes through to delta_maps.py.
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
  PYTHONPATH=$CL:$CAMPAIGN/jobs OMP_NUM_THREADS=1 $CL/.venv/bin/python delta_maps.py $m \
    $CAMPAIGN/runs/delta/${m}${SUFFIX:-_p5}.csv
done
