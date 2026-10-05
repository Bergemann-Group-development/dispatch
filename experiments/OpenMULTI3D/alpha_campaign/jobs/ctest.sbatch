#!/bin/bash -l
# The crop agreement test (make_nml.py crop_test): ctest.sbatch 1 runs the 6x6 3D NLTE
# runs (Fe cropped and uncropped, Mg) and the plain LTE run side by side; ctest.sbatch 2
# links and maps the departures and runs the three stage-2 SEDs. NLTE_BIN is the NLTE
# build (needs save_every, 1631491e2), SED_BIN the production SED build, HEAD_BIN the
# build that sed_full_head checks. (common.sh has already set BIN to bin/current.)
#   sbatch --nice=100 -w mpia-slurm-node88 -c 20 --mem=48G jobs/ctest.sbatch 1
#   sbatch --nice=100 -w mpia-slurm-node88 -c 18 --mem=48G -d afterok:JOB jobs/ctest.sbatch 2
#SBATCH -J ctest
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/ctest-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
NLTE_BIN=${NLTE_BIN:-$CAMPAIGN/bin/4890f9fde}
SED_BIN=${SED_BIN:-$CAMPAIGN/bin/23675093a}
HEAD_BIN=${HEAD_BIN:-$CAMPAIGN/bin/4890f9fde}
T=ctest/t60g40m20_s00091
export TIMINGS=$CAMPAIGN/logs/timings_ctest_$SLURM_JOB_ID.txt
echo "host $(hostname -s), stage $1, $SLURM_CPUS_PER_TASK cores"
cd $WORK
if [ "$1" = 1 ]; then
  BIN=$NLTE_BIN
  run $T/fe_full.nml 8 & run $T/fe_crop.nml 6 & run $T/mg_full.nml 4 & run $T/lte.nml 2 &
  wait
else
  python3 $CAMPAIGN/jobs/make_nml.py crop_test_dep || exit 1
  R=$CAMPAIGN/runs/$T
  [ -f $R/dep_crop/save/pop_Fe_meta.txt ] || \
    PYTHONPATH=$CL OMP_NUM_THREADS=1 $CL/.venv/bin/python $CAMPAIGN/jobs/map_departures.py \
      $R/fe_crop $R/mg_full $R/dep_crop Fe || exit 1
  BIN=$SED_BIN run $T/sed_full.nml 6 & BIN=$SED_BIN run $T/sed_crop.nml 6 &
  BIN=$HEAD_BIN run $T/sed_full_head.nml 6 &
  wait
fi
cat $TIMINGS
