# Default for gfortran
OPT = -O1
HW =
W = -Wline-truncation -Wno-surprising -w -fno-range-check
DEBUG = -g -fbacktrace 

# Magic strings for dependencies
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC  = mpifort -I${BUILD} -J${BUILD}
else
  FC  = gfortran -I${BUILD} -J${BUILD}
endif
F90 = $(FC) -x f95-cpp-input -ffree-line-length-0 $(W)
LD  = $(FC)

# Set openmp compiler flags if needed
ifdef OMP
ifeq ($(strip $(OMP)),yes)
  PAR = -fopenmp
endif
endif

FFLAGS = $(OPT) $(PAR) $(HW) $(DEBUG)
LDFLAGS = $(OPT) $(PAR) $(DEBUG)

