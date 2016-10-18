# Default for cray fortran.
# Note cray does not support auto-generation of dependencies (AUTO=no).
OPT = -O1
DEBUG =
HW =
AUTO = no

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  F90 = ftn -eZ
  LD  = ftn
else
  F90 = ftn -target-network=none -eZ
  LD  = ftn -target-network=none
endif

# Set openmp compiler flags if needed
ifeq ($(strip $(OMP)),yes)
PAR = -homp
endif

# Pick up sub-options in the right order
FFLAGS  = $(OPT) $(HW) $(PAR) $(DEBUG)
LDFLAGS = $(OPT) $(HW) $(PAR) $(DEBUG)
