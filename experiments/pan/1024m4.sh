#!/bin/bash
#SBATCH --job-name=pan_1024m4
#SBATCH --ntasks=32 --ntasks-per-node=1 --ntasks-per-core=1 --cpus-per-task=48
##SBATCH --time=12:00:00

export S=${HOME}/codes/dispatch2/utilities/scripts

$S/jobcontrol.csh $S/slurm/intel.sh 1024m4 dispatch_stagger2e.x
