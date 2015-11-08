VPATH := $(VPATH):$(TOP)/$(MPI)

# Local object files to include
#OBJ += \
#omp_mod.o \
#mpi_basic.o \
#mpi_all.o \
#mpi_coords.o \
#mpi_io.o \
#mpi_reduce.o \
#mpi_send.o \
#mpi_rma.o

mpi_basic.o:  omp_mod.o
mpi_all.o:    mpi_basic.o
mpi_coords.o: mpi_basic.o
mpi_io.o:     mpi_basic.o omp_mod.o
mpi_reduce.o: mpi_basic.o
mpi_rma.o:    mpi_basic.o mpi_io.o
mpi_send.o:   mpi_basic.o
mpi_mesg.o:   mpi_basic.o

