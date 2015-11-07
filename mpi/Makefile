VPATH := $(VPATH):$(TOP)/$(MPI)

# Local object files to include
OBJ += \
omp_mod.o \
mpi_base.o

#mpi_all.o \
#mpi_coords.o \
#mpi_io.o \
#mpi_reduce.o \
#mpi_send.o \
#mpi_rma.o

mpi_base.o:   omp_mod.o
mpi_all.o:    mpi_base.o
mpi_coords.o: mpi_base.o
mpi_io.o:     mpi_base.o omp_mod.o
mpi_reduce.o: mpi_base.o
mpi_rma.o:    mpi_base.o mpi_io.o
mpi_send.o:   mpi_base.o
mpi_mesg.o:   mpi_base.o

