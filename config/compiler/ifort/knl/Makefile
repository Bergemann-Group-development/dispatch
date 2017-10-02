HW    = -xMIC-AVX512
PAR   = -qopenmp -qopt-threads-per-core=1

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC  = mpiifort
else
  FC  = ifort
endif
