# Default for ifort. Beware that versions prior to version 12.3 may crash
OPT = -O1
DEBUG = -g -traceback -fpe:0
#HW = -axAVX -xSSE4.1
HW =

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  F90 = mpifort -cpp -assume byterecl -diag-disable remark
  LD  = mpifort
else
  F90 = ifort -cpp -assume byterecl -diag-disable remark
  LD  = ifort
endif

# Set openmp compiler flags if needed
ifeq ($(strip $(OMP)),yes)
PAR = -qopenmp
endif

# Pick up sub-options in the right order
FFLAGS  = $(DEBUG) $(OPT) $(HW) $(PAR)
LDFLAGS = $(DEBUG) $(OPT) $(HW) $(PAR)
