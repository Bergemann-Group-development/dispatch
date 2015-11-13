VPATH := $(VPATH):$(TOP)/$(MPI)

# Dispatch object files to include
ifeq ($(EXEC),dispatch.x)
  OBJ += \
  omp_mod.o \
  mpi_timer.o \
  mpi_basic.o \
  mpi_coords.o
  #mpi_all.o \
  #mpi_reduce.o \
  #mpi_send.o \
  #mpi_rma.o
endif

mpi_basic.o:  omp_mod.o
mpi_all.o:    mpi_basic.o
mpi_coords.o: mpi_basic.o
mpi_par_io.o: mpi_basic.o omp_mod.o mpi_timer.o
mpi_reduce.o: mpi_basic.o
mpi_rma.o:    mpi_basic.o
mpi_send.o:   mpi_basic.o
mpi_mesg.o:   mpi_basic.o

