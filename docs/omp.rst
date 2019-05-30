OpenMP
========

OpenMPI is fundamentatl in DISPATCH, so to take advantage of
task parallelization on a single node, set the OMP_NUM_THREADS
environemental variable before starting. For example:::

  export OMP_NUM_THREADS=20
  ./dispatch.x

To see how many threads your compute node supports::

  grep -c processor /proc/cpuinfo

That file contains other detailed information about the cores,
which might be relevant when choosing compiler options ::

  more /proc/cpuinfo

