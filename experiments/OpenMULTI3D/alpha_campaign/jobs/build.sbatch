#!/bin/bash -l
# Build dispatch.x, libhop.so and hop_prune.x with ifort on a compute node, then
# freeze copies under $CAMPAIGN/bin/<short hash>/ so later rebuilds cannot change the
# binaries a queued campaign job runs.  CLONE defaults to the campaign's own checkout
# $CAMPAIGN/src (fixes are applied there with git am); bin/current moves only with
# SET_CURRENT=1.  Size -c and -j to the free cores, e.g. sbatch -c 1 --mem=4G for an
# incremental build on a nearly full node.
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
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
CLONE=${CLONE:-$CAMPAIGN/src/experiments/OpenMULTI3D}
cd "$CLONE"
echo "host $(hostname)"; lscpu | grep -E 'Model name|^CPU\(s\)|Socket'
lscpu | grep -q avx512f && echo "avx512: yes" || echo "avx512: no"
which ifort; ifort --version 2>&1 | head -1 || true
rev=$(git rev-parse --short=9 HEAD)
echo "revision $rev"; git status --short | head
/usr/bin/time -v make -j$SLURM_CPUS_PER_TASK COMPILER=ifort > "$CAMPAIGN/logs/make-$SLURM_JOB_ID.log" 2>&1; rc=$?; tail -40 "$CAMPAIGN/logs/make-$SLURM_JOB_ID.log"; [ $rc -eq 0 ] || { echo "BUILD FAILED rc=$rc"; exit 1; }
ls -l dispatch.x HOP/standalone/build/libhop.so HOP/standalone/build/hop_prune.x
mkdir -p "$CAMPAIGN/bin/$rev"
cp -p dispatch.x HOP/standalone/build/libhop.so HOP/standalone/build/hop_prune.x "$CAMPAIGN/bin/$rev/"
# REVISION names the commit and tree of the sources, which jobs/common.sh logs. The start
# banner of dispatch.x prints the same hash (the Makefile rebuilds it when HEAD moves).
echo "$(git rev-parse HEAD) tree $(git rev-parse 'HEAD^{tree}')$([ -n "$(git status --porcelain --untracked-files=no)" ] && echo ' (modified)')" > "$CAMPAIGN/bin/$rev/REVISION"
[ "${SET_CURRENT:-0}" = 1 ] && ln -sfn "$rev" "$CAMPAIGN/bin/current"
echo "frozen in $CAMPAIGN/bin/$rev"
