#!/bin/bash -l
# queue.sbatch LIST [THREADS]: work through a task list (relative to work/), one task per
# line, until every line is finished or claimed.
#
# The job runs cpus-per-task/THREADS workers (THREADS defaults to all the job's cores),
# so a job sized to the free cores of a partial node uses exactly those. Any number of
# jobs and array tasks can serve one list: a task is claimed by an atomic mkdir under
# CLAIMS (default logs/claims_<list>), and finished or claimed tasks are skipped. A claim
# whose job has ended without finishing the task is taken over by the next job, unless
# the task failed; RETRY_FAILED=1 retries those too. Each job reads a private copy of
# the list made at its start, so the list may be edited or regenerated while jobs run.
#
# By default a line is a namelist, run with dispatch.x and finished once its log says
# "OpenMULTI3D terminated". STEPS=name sources jobs/steps/name.sh, which may redefine
#   task_key LINE...   claim name of a line
#   task_done LINE...  true if the line needs no work
#   task LINE...       do the work
#   prepare            run once per submission, before the list is copied
#   sbatch -c 4 --mem=30G -w mpia-slurm-node83 jobs/queue.sbatch x3d/tasks.txt
#   sbatch --array=0-3 --exclusive -c 32 --mem=110G -t 6:00:00 jobs/queue.sbatch b1/tasks.txt 1
#   sbatch -c 1 --mem=4G --export=ALL,STEPS=dep4,BIN=$C/bin/23675093a jobs/queue.sbatch b4/tasks_dep.txt
# The OpenMP runtime comes from lib/libiomp5.so (copied from the login node's oneAPI
# 2024.0), because /opt/intel is a stale mount on node83.
#SBATCH -J alpha_queue
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --nodelist=mpia-slurm-node[81-90]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/logs/queue-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
export LD_LIBRARY_PATH=$CAMPAIGN/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
LIST=$1
NT=${2:-$SLURM_CPUS_PER_TASK}
P=$(( SLURM_CPUS_PER_TASK / NT )); [ $P -ge 1 ] || P=1
tag=$(echo "${LIST%.txt}" | tr / _)
ME=${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID}
ME=${ME:-$SLURM_JOB_ID}
CLAIMS=${CLAIMS:-$CAMPAIGN/logs/claims_$tag}
export TIMINGS=$CAMPAIGN/logs/timings_${tag}_$ME.txt
mkdir -p "$CLAIMS"
cd "$WORK"

task_key () { local n=${1%.nml}; echo "${n//\//_}"; }
task_done () { grep -qs 'OpenMULTI3D terminated' "$CAMPAIGN/logs/$(task_key "$1").log"; }
task () { run "$1" "$NT"; grep "^${1%.nml} " "$TIMINGS" | tail -1; task_done "$1"; }
prepare () { :; }
[ -n "${STEPS:-}" ] && source "$CAMPAIGN/jobs/steps/$STEPS.sh"

# claim KEY: true if this job now owns the task. A claim records its owner; one left by
# an ended job is moved aside (atomically, so only one claimant wins) and made anew.
claim () {
  local d=$CLAIMS/$1 owner live
  if mkdir "$d" 2>/dev/null; then echo "$ME" > "$d/owner"; return 0; fi
  [ -f "$d/failed" ] && [ "${RETRY_FAILED:-0}" != 1 ] && return 1
  owner=$(cat "$d/owner" 2>/dev/null) || return 1   # claims without owner are kept
  live=$(squeue -h -u "$USER" -o %i) || return 1
  grep -qx "$owner" <<< "$live" && return 1
  mv "$d" "$d.stale.$ME.$BASHPID" 2>/dev/null || return 1
  mkdir "$d" 2>/dev/null || return 1
  echo "$ME" > "$d/owner"
  echo "took over $1 from ended job $owner"
}

# one submission (array job or single job) runs prepare once; the others wait for it
if mkdir "$CLAIMS/.prepare_${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}" 2>/dev/null; then
  prepare || { echo "prepare failed"; exit 1; }
  touch "$CLAIMS/.prepare_${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}/done"
else
  until [ -f "$CLAIMS/.prepare_${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}/done" ]; do sleep 10; done
fi
SNAP=$(mktemp) && cp "$WORK/$LIST" "$SNAP" || exit 1

worker () {
  local line key
  while read -r line; do
    [ -z "$line" ] || [ "${line:0:1}" = '#' ] && continue
    task_done $line && continue
    key=$(task_key $line)
    claim "$key" || continue
    task $line || { echo "$key failed"; touch "$CLAIMS/$key/failed"; }
  done < "$SNAP"
}

echo "host $(hostname -s), $P workers x $NT threads, list $LIST, BIN=$BIN${STEPS:+, steps $STEPS}"
for i in $(seq $P); do worker & done
wait
rm -f "$SNAP"
[ -f "$TIMINGS" ] && echo "done: $(wc -l < "$TIMINGS") runs, $(awk '$5!=0' "$TIMINGS" | wc -l) nonzero exits"
