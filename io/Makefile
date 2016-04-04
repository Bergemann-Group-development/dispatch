# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Dependencies
io_unit_mod.o:
io_mod.o:               omp_mod.o omp_timer_mod.o mpi_mod.o io_unit_mod.o
trace_mod.o:            io_mod.o
buffered_output_mod.o:  trace_mod.o
