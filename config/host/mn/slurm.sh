#!/bin/bash
#SBATCH --qos=class_a --job-name=jobname
#SBATCH --ntasks=256 --ntasks-per-node=2 --ntasks-per-core=1 --cpus-per-task=24
## Max number of nodes (= ntasksi*cpus-per-task/48) is 200 in class_a

# Path to scripts
export S=$HOME/codes/dispatch2/utilities/scripts

# Example, with input file 2048m4.nml and output to data/2048m4
$S/jobcontrol.csh $S/slurm/intel.sh 2048m4 dispatch.x
