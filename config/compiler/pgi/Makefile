# Default for PGI
OPT = -O2
HW =
W = -Mpreprocess
DEBUG =

# Magic strings for dependencies
DEPFLAGS =
POSTCOMPILE =

# Depending on MPI use MPI-option or not
FC  = pgfortran
ifeq ($(strip $(MPI)),mpi)
  MPIOPT = -Mmpi=msmpi
else
  MPIOPT =
endif
F90 = $(FC) $(MPIOPT) $(W)
LD  = $(FC) $(MPIOPT)

# Set openmp compiler flags if needed
ifdef OMP
ifeq ($(strip $(OMP)),yes)
PAR =
endif
endif

FFLAGS = $(OPT) $(PAR) $(HW) $(DEBUG)
LDFLAGS = $(OPT) $(PAR) $(DEBUG)
