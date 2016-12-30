#!/bin/bash
#SBATCH --partition=astro_devel
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2 --cpus-per-task=10
##SBATCH --ntasks-per-node=2 --cpus-per-task=10 --threads-per-core=2
##SBATCH --ntasks-per-node=1 --cpus-per-task=20
##SBATCH --ntasks-per-core=2 --cpus-per-task=20
#SBATCH --job-name=pa

cd $SLURM_SUBMIT_DIR
/bin/rm -f data/devel/*.log

# time stamp files
run="data"
time=`date +%y%m%d_%H%M`
log=$run/log_$time
x=./dispatch.x
ln -sf $log log

# redirect stdout and stderr
exec >$log 2>&1

# number of threads and OpenMP affinity
export KMP_AFFINITY="granularity=core,scatter"
export OMP_NUM_THREADS=${1:-$SLURM_CPUS_PER_TASK}
echo "SLURM_JOB_ID    = $SLURM_JOB_ID"
echo "KMP_AFFINITY    = $KMP_AFFINITY"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"

cat >devel.nml <<EOF
&experiment_params /
&mpi_mesg_params   max_sent=500 max_recv=100 min_nq=10 /
&timer_params      sec_per_report=10 /
&io_params         verbose=1 do_trace=f do_debug=f do_output=f omp_trace=t /
&box_params        size=4,2,2 dims=40,20,20 mpi_dims=4,2,2 origin=0,0,0 /
&patch_params      nt=5 nw=1 courant=0.2 grace=0.3 n=3*32 /
&out_params        end_time=0.02 out_time=1 print_time=0 guard_zones=f time_derivs=f /
&force_params      k=1,1,2 a0=3,3.2,3.1 t_turn=0.30 type='single_solenoidal' /
&initial_params    k=3,3,3 a0=0.0,0.0,0.0 u0=0,0,0 b0=0,0,0.0 type='single_compressive' /
&ramses_params     detailed_timer=f courant_factor=0.8 gamma=1.00001 slope_type=2 smallr=1e-4 smallc=5e-1 riemann='hllc' do_isothermal=t /
EOF

export OMP_NUM_THREADS=10
srun --label --distribution=block:block ./dispatch-mpi-1-3-2.x devel.nml
