# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Dependencies
io_unit_mod.o:
os_mod.o:               io_unit_mod.o mpi_mod.o omp_mod.o
io_mod.o:               omp_mod.o omp_timer_mod.o mpi_mod.o io_unit_mod.o timer_mod.o os_mod.o
trace_mod.o:            io_mod.o timer_mod.o
timer_mod.o:            omp_mod.o omp_timer_mod.o omp_lock_mod.o io_unit_mod.o mpi_mod.o
legacy_io_mod.o:        io_mod.o trace_mod.o patch_mod.o kinds_mod.o
buffered_io_mod.o:      io_mod.o trace_mod.o patch_mod.o
direct_io_mod.o:        io_mod.o trace_mod.o patch_mod.o
parallel_io_mod.o:      io_mod.o trace_mod.o patch_mod.o mpi_io_mod.o central_mod.o
snapshot_list_mod.o:    io_mod.o trace_mod.o patch_mod.o 
data_io_mod.o:          io_mod.o trace_mod.o patch_mod.o direct_io_mod.o parallel_io_mod.o buffered_io_mod.o legacy_io_mod.o \
                        snapshot_list_mod.o
