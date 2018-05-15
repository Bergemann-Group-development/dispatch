VPATH := $(VPATH):$(TOP)/par/mpi

mpi_mod.o:         omp_mod.o io_unit_mod.o omp_timer_mod.o
mpi_coords_mod.o:  mpi_mod.o io_mod.o
mpi_mesg_mod.o:    mpi_mod.o omp_mod.o omp_timer_mod.o omp_lock_mod.o trace_mod.o io_mod.o
mpi_io_file_mod.o: mpi_mod.o io_mod.o trace_mod.o
mpi_io_mod.o:      mpi_mod.o io_mod.o trace_mod.o mpi_io_file_mod.o
