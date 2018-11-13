#!/bin/bash
#SBATCH --job-name=pan_2048m4
#SBATCH --ntasks=256 --ntasks-per-node=2 --ntasks-per-core=1 --cpus-per-task=24

export S=${HOME}/codes/dispatch2/utilities/scripts

$S/jobcontrol.csh $S/slurm/intel.sh 2048m4 dispatch_stagger2e.x
