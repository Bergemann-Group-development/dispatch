VPATH := $(VPATH):$(TOP)/mpi

omp_timer_mod.o:  omp_mod.o
omp_lock_mod.o:   omp_mod.o
mpi_mod.o:        omp_mod.o io_mod.o
mpi_all_mod.o:    mpi_mod.o
mpi_coords_mod.o: mpi_mod.o io_mod.o
mpi_par_io_mod.o: mpi_mod.o mpi_coords_mod.o mpi_timer.o io_mod.o
mpi_reduce_mod.o: mpi_mod.o io_mod.o
mpi_rma_mod.o:    mpi_mod.o
mpi_send_mod.o:   mpi_mod.o
mpi_mesg_mod.o:   mpi_mod.o

