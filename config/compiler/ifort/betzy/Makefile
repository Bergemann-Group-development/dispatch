OPT     = -O3
HW      = -xavx2 -qopt-prefetch -align array64byte -qopt-assume-safe-padding -fma -fp-model source -qopenmp
DEBUG   = -traceback
LDFLAGS = -qopenmp

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC  = mpif90
else
  FC  = ifort
endif
