#!/bin/bash -l
# Build dispatch.x, libhop.so and hop_prune.x with ifort on a compute node, then
# freeze copies under $CAMPAIGN/bin/<short hash>/ so later rebuilds of the clone
# cannot change the binaries a queued campaign job runs.
#SBATCH -J alpha_build
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/build-%j.out
set -uo pipefail
CLONE=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/dispatch/experiments/OpenMULTI3D
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
cd "$CLONE"
echo "host $(hostname)"; lscpu | grep -E 'Model name|^CPU\(s\)|Socket'
lscpu | grep -q avx512f && echo "avx512: yes" || echo "avx512: no"
which ifort; ifort --version 2>&1 | head -1 || true
rev=$(git rev-parse --short=9 HEAD)
echo "revision $rev"; git status --short | head
/usr/bin/time -v make -j16 COMPILER=ifort > "$CAMPAIGN/logs/make-$SLURM_JOB_ID.log" 2>&1; rc=$?; tail -40 "$CAMPAIGN/logs/make-$SLURM_JOB_ID.log"; [ $rc -eq 0 ] || { echo "BUILD FAILED rc=$rc"; exit 1; }
ls -l dispatch.x HOP/standalone/build/libhop.so HOP/standalone/build/hop_prune.x
mkdir -p "$CAMPAIGN/bin/$rev"
cp -p dispatch.x HOP/standalone/build/libhop.so HOP/standalone/build/hop_prune.x "$CAMPAIGN/bin/$rev/"
ln -sfn "$rev" "$CAMPAIGN/bin/current"
echo "frozen in $CAMPAIGN/bin/$rev"
