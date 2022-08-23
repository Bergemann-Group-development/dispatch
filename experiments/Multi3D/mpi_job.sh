#!/bin/bash -l
# --- Log file and error file names ---
regex="datadir='\K(.*?)(?=')"
datadir=$(grep -oP $regex $1)
regex="[\w\/\.]+(?=[\.])"
rundir=$(grep -oP $regex <<< $1)
if [ $rundir != "input" ]
then
datadir=$datadir$rundir'/'
fi
echo $datadir
mkdir -p $datadir
# -------------------------------------
sbatch -o $datadir'slurm.log' -e $datadir'slurm.err' 2nodes.slurm $1