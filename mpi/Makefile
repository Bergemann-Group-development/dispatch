VPATH := $(VPATH):$(TOP)/mpi

mpi_mod.o:         omp_mod.o io_unit_mod.o omp_timer_mod.o
mpi_coords_mod.o:  mpi_mod.o io_mod.o
mpi_mesg_mod.o:    mpi_mod.o io_mod.o io_unit_mod.o omp_mod.o omp_timer_mod.o omp_lock_mod.o trace_mod.o
mpi_file_mod.o:    mpi_mod.o io_mod.o trace_mod.o
mpi_io_mod.o:      mpi_mod.o io_mod.o io_unit_mod.o trace_mod.o mpi_file_mod.o omp_mod.o \
                   omp_timer_mod.o dll_mod.o
mpi_comm_mod.o:    mpi_mod.o
mpi_rma_mod.o:     mpi_mod.o
mpi_buffer_mod.o:  trace_mod.o
mpi_counter_mod.o: mpi_mod.o io_unit_mod.o omp_timer_mod.o trace_mod.o
mpi_global_mod.o:  mpi_mod.o trace_mod.o
