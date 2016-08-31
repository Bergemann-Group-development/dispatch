#PBS -W group_list=g26202
##PBS -W group_list=g26161       # Bob
#PBS -N bro
#PBS -q devel
#PBS -S /bin/csh      
#PBS -l select=1:ncpus=16:mpiprocs=2:ompthreads=8:model=bro
#PBS -l walltime=0:05:00
#PBS -j oe
#PBS -m e

#module load comp-intel/11.1.072 mpi-sgi/mpt.2.04.10789

csh test.sh
