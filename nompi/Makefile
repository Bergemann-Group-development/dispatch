VPATH := $(VPATH):$(TOP)/$(MPI)

# Dispatch object files to include
ifeq ($(EXEC),dispatch.x)
  OBJ += \
  omp_mod.o \
  mpi_basic.o \
  mpi_reduce.o \
  mpi_io.o
endif

mpi_basic.o:  omp_mod.o
mpi_reduce.o: mpi_basic.o
mpi_io.o:     mpi_basic.o

