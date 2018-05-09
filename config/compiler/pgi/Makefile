# Default for PGI
OPT = -O2
HW =
W =
DEBUG =

# Magic strings for dependencies
DEPFLAGS =
POSTCOMPILE =

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
FC  = pgfortran -Mpreprocess -Mmpi=msmpi $(W)
else
FC  = pgfortran -Mpreprocess $(W)
endif
F90 = $(FC) $(W)
LD  = $(FC)

# Set openmp compiler flags if needed
ifdef OMP
ifeq ($(strip $(OMP)),yes)
PAR =
endif
endif

FFLAGS = $(OPT) $(PAR) $(HW) $(DEBUG)
LDFLAGS = $(OPT) $(PAR) $(DEBUG)
