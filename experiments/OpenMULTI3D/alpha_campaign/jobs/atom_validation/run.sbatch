#!/bin/bash -l
# Atom validation: 1D NLTE runs of a task list (make_runs.py), PAR at a time with NT
# threads each, in $C/runs/atom_validation. Submit with --export=ALL,LIST=tasks_<batch>.txt
# and BIN set to the frozen build.
#SBATCH -J atom_valid
#SBATCH --partition=four-wks
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/runs/atom_validation/logs/job-%j.out
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
RUNS=$CAMPAIGN/runs/atom_validation
LIST=${LIST:?task list}
PAR=${PAR:-16}
NT=${NT:-2}
cd "$RUNS"
mkdir -p logs
# short relative paths the namelists use
ln -sfn "$WORK/input" input
ln -sfn "$WORK/absmet" absmet
ln -sfn "$CAMPAIGN/runs/atom_audit/atoms" atoms
ln -sfn "$CAMPAIGN/scratch/atom_audit/public" public
bin_info
ulimit -s unlimited
one () {
  local nml=$1 name; name=${nml#nml/}; name=${name%.nml}; name=${name//\//_}
  OMP_NUM_THREADS=$NT /usr/bin/time -f "$name %e %M %x $(hostname -s)" -a -o "$RUNS/logs/timings.txt" \
    "$BIN/dispatch.x" "$nml" > "$RUNS/logs/$name.log" 2>&1
}
export -f one; export RUNS NT BIN KMP_STACKSIZE
xargs -P "$PAR" -I{} bash -c 'one {}' < "$LIST"
echo done
