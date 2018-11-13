#!/bin/bash
#SBATCH --qos=debug --job-name=debug
#SBATCH --ntasks=8 --ntasks-per-node=1 --ntasks-per-core=1 --cpus-per-task=48

export S=${HOME}/codes/dispatch2/utilities/scripts
$S/jobcontrol.csh $S/slurm/intel.sh input dispatch.x
