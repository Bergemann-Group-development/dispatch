Running under MPI
=================

Use the relevant cluster documentation to find out how to run
under MPI.  To use SLURM to run a job with 64 MPI processes,
with two processes per node and 10 cores per process, you
may need something similar to this (which would work on HPC.KU:DK)::

  #!/bin/bash
  #SBATCH --ntasks=64 --ntasks-per-node=2 --cpus-per-task=10
  #SBATCH --mem-per-cpu=3g --exclusive

  export OMP_NUM_THREADS=10
  ./dispatch.x

Hyper-threading
---------------

If the system supports hyper-threading; e.g., 2 threads per core, with
10 cores per process:::

  #!/bin/bash
  #SBATCH --ntasks=64 --ntasks-per-node=2 --cpus-per-task=10 --ntasks_per_core=2
  #SBATCH --mem-per-cpu=3g --exclusive
  
  export OMP_NUM_THREADS=20
  ./dispatch.x
