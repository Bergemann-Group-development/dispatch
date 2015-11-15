VPATH := $(VPATH):$(TOP)/$(MPI)

# Dispatch object files to include
ifeq ($(EXEC),dispatch.x)
  OBJ += \
  omp_mod.o \
  mpi_timer.o \
  mpi_mod.o \
  mpi_coord_mod.o
  #mpi_all.o \
  #mpi_reduce.o \
  #mpi_send.o \
  #mpi_rma.o
endif

mpi_mod.o:        omp_mod.o
mpi_all_mod.o:    mpi_mod.o
mpi_coord_mod.o:  mpi_mod.o
mpi_par_io_mod.o: mpi_mod.o mpi_coord_mod.o mpi_timer.o
mpi_reduce_mod.o: mpi_mod.o
mpi_rma_mod.o:    mpi_mod.o
mpi_send_mod.o:   mpi_mod.o
mpi_mesg_mod.o:   mpi_mod.o

