# Default for ifort. Beware that versions prior to version 12.3 may crash
OPT = -O1
DEBUG = -g -traceback -fpe:0
#HW = -axAVX -xSSE4.1
HW =
PREC = 4

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  F90 = mpiifort -cpp -assume byterecl -diag-disable remark -DPREC=$(PREC)
  LD  = mpiifort
else
  F90 = ifort -cpp -assume byterecl -diag-disable remark -DPREC=$(PREC)
  LD  = ifort
endif

# Set openmp compiler flags if needed
ifeq ($(strip $(OMP)),yes)
PAR = -qopenmp
endif

# Pick up sub-options in the right order
FFLAGS  = $(DEBUG) $(OPT) $(HW) $(PAR)
LDFLAGS = $(DEBUG) $(OPT) $(HW) $(PAR)
