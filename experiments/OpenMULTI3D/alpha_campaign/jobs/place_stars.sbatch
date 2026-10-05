#!/bin/bash -l
# Crude 1D LTE labels of the 4MOST HRS spectra (place_stars.py): unpack obs/*.zip, put
# the batch-1a grid on the observed wavelengths, fit every spectrum.  Spare cores only.
#   sbatch --nodelist=mpia-slurm-node83 jobs/place_stars.sbatch
#SBATCH -J place_stars
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/place-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
export PYTHONPATH=$CLONE OMP_NUM_THREADS=1 HOP_LIB=$BIN/libhop.so
OUT=$CAMPAIGN/runs/place
mkdir -p "$OUT"
cd "$CAMPAIGN/obs"
for z in 0*.zip; do [ -d "${z%.zip}" ] || unzip -q "$z" -x '__MACOSX/*'; done
PY="$CLONE/.venv/bin/python $CAMPAIGN/jobs/place_stars.py"
[ -f "$OUT/grid_R20000.npz" ] || $PY grid "$OUT/grid_R20000.npz"
$PY fit "$OUT/grid_R20000.npz" "$OUT/labels.csv"
