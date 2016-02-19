VPATH := $(VPATH):$(TOP)/$(MPI)

mpi_basic.o:  omp_mod.o
mpi_reduce.o: mpi_basic.o
mpi_io.o:     mpi_basic.o

