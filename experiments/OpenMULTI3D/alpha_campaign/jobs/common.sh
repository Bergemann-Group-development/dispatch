# Shared settings of the alpha-campaign job scripts (sourced, not run).
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
WORK=$CAMPAIGN/work
BIN=${BIN:-$CAMPAIGN/bin/current}    # a job may export another frozen build
export KMP_STACKSIZE=256M
unset OMP_PROC_BIND OMP_PLACES I_MPI_PMI_LIBRARY

# run NML THREADS: run dispatch.x on WORK/NML with THREADS OpenMP threads; append
# "name threads wall_s maxrss_kb exit host" to $TIMINGS and keep the log beside it.
run () {
  local nml=$1 nt=$2 name=${1%.nml}
  local log=$CAMPAIGN/logs/${name//\//_}.log
  ( cd "$WORK" && OMP_NUM_THREADS=$nt /usr/bin/time -f "%e %M %x" -o "$log.time" \
      "$BIN/dispatch.x" "$nml" > "$log" 2>&1 )
  read -r wall rss rc < "$log.time"
  echo "$name $nt $wall $rss $rc $(hostname -s)" >> "$TIMINGS"
}

# variant NML TAG: copy WORK/NML to a sibling NAME_TAG.nml (its own output folder)
variant () { cp "$WORK/$1" "$WORK/${1%.nml}_$2.nml"; echo "${1%.nml}_$2.nml"; }

# The campaign may only use nodes 81-90 (booked); a job that lands elsewhere stops here.
if [ -n "$SLURM_JOB_ID" ]; then
  case $(hostname -s) in
    mpia-slurm-node8[1-9]|mpia-slurm-node90) ;;
    *) echo "$(hostname -s) is outside nodes 81-90; stopping"; exit 1 ;;
  esac
fi
