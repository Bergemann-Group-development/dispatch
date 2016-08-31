#PBS -N ivy
#PBS -q devel
#PBS -S /bin/csh      
#PBS -l select=1:ncpus=16:mpiprocs=2:ompthreads=8:model=ivy
#PBS -l walltime=0:05:00
#PBS -j oe
#PBS -m e

#module load comp-intel/11.1.072 mpi-sgi/mpt.2.04.10789

#setenv MPI_BUFS_PER_HOST 1024
#setenv MPI_BUFS_PER_PROC 256

csh test.sh
