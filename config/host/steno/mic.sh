export I_MPI_ENV_PREFIX_LIST=knc:MIC
export MIC_OMP_NUM_THREADS=$2
export OMP_NUM_THREADS=5
export KMP_AFFINITY=granularity=fine,scatter,0,1
mpiexec.hydra -machine mics -ppn 1 -l -n $1 $3 $4
