# Set standard values to use gfortran


# Default for gfortran
OPT = -O2
HW =
PREC = 4
W = -Wline-truncation -Wno-surprising -w
D = -DPREC=$(PREC)

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
F90 = mpifort -x f95-cpp-input -ffree-line-length-0 $(W) $(D)
LD  = mpifort
else
F90 = gfortran -x f95-cpp-input -ffree-line-length-0 $(W) $(D)
LD  = gfortran
endif

# Set openmp compiler flags if needed
ifdef OMP
ifeq ($(strip $(OMP)),yes)
PAR = -fopenmp
endif
endif

FFLAGS = $(OPT) $(PAR) $(HW) $(DEBUG)
LDFLAGS = $(OPT) $(PAR) $(DEBUG)
