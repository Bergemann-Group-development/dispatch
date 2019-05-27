# Default for ifort. Beware that versions prior to version 12.3 may crash
OPT = -O1
DEBUG = -g -traceback
#HW = -axAVX -xSSE4.1
HW =

# Magic strings for dependencies
DEPFLAGS = -gen-depformat=make -gen-dep=$(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC = mpiifort -module $(BUILD)
else
  FC = ifort -module $(BUILD)
endif
F90 = $(FC) -cpp -assume byterecl -diag-disable remark
LD  = $(FC)

# Set openmp compiler flags if needed
ifeq ($(strip $(OMP)),yes)
  PAR = -qopenmp -parallel-source-info=2
endif

# Pick up sub-options in the right order
FFLAGS  = $(DEBUG) $(OPT) $(HW) $(PAR)
LDFLAGS = $(DEBUG) $(OPT) $(HW) $(PAR)
