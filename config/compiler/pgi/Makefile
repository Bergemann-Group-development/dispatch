# Default for PGI
OPT = -O2
HW =
W = -Mpreprocess -fast -Mipa=fast,inline -llapack -lblas -tp=haswell-64 -traceback -g -Mcpp=mmd
#DEBUG = -g -traceback -fpe:0

# Magic strings for dependencies
DEPFLAGS =
POSTCOMPILE =

# Depending on MPI use MPI-option or not
FC  = mpifort
ifeq ($(strip $(MPI)),mpi)
 MPIOPT = 
else
 MPIOPT = 
endif
F90 = $(FC) $(MPIOPT) $(W)
LD  = $(FC) $(MPIOPT) $(W)

# Set openmp compiler flags if needed
ifdef OMP
ifeq ($(strip $(OMP)),yes)
PAR = -omp
endif
endif

FFLAGS = $(OPT) $(PAR) $(HW) $(DEBUG)
LDFLAGS = $(OPT) $(PAR) $(DEBUG)
