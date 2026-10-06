#!/bin/bash -l
# avg3d_queue.sbatch AVG [RECIPE]: the <3D> average of every model in work/AVG/models.txt
# (make_nml.py avg3d: 'name [Fe/H] [alpha/Fe] snapshot ...') into work/AVG/<model>.atmos
# with avg3d.py (RECIPE t4 or geo, default t4).  Claim-based (logs/claims_<AVG>), so
# several runners can share a list; runs under Slurm and on the vis nodes:
#   sbatch -w NODE -c 4 --mem=16G --export=ALL,BIN=$C/bin/<hash> jobs/avg3d_queue.sbatch avg3d8
#   ssh vis2 "cd $C && BIN=$C/bin/<hash> nohup jobs/avg3d_queue.sbatch avg3d8 > logs/avg3d8_vis2.out 2>&1 &"
# tumult comes from PY (default the campaign checkout $C/src, which must hold a commit
# with tumult.avg3d's top-of-column options) and libhop.so from BIN.  About 30 s and
# 6 GB per model.
#SBATCH -J alpha_avg3d
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/avg3d-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CL=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
PY=${PY:-$CAMPAIGN/src/experiments/OpenMULTI3D}
AVG=$1
RECIPE=${2:-t4}
CLAIMS=${CLAIMS:-$CAMPAIGN/logs/claims_$AVG}
mkdir -p "$CLAIMS"
ulimit -s unlimited
set -o pipefail
cd "$WORK"
while read -r model feh alpha snaps; do
  [ -z "$model" ] && continue
  [ -s "$AVG/$model.atmos" ] && continue
  mkdir "$CLAIMS/$model" 2>/dev/null || continue
  echo "$(hostname -s) ${SLURM_JOB_ID:-$$}" > "$CLAIMS/$model/owner"
  HOP_LIB=$BIN/libhop.so HOP_HEADER=$PY/HOP/standalone/libhop.h PYTHONPATH=$PY \
    OMP_NUM_THREADS=${OMP_NUM_THREADS_PY:-4} $CL/.venv/bin/python $CAMPAIGN/jobs/avg3d.py \
    "$AVG/$model.atmos" "$feh" "$alpha" $snaps --recipe "$RECIPE" | cut -c1-160 \
    || echo "$model: average failed"
done < "$WORK/$AVG/models.txt"
