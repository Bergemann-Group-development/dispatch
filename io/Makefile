# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Dependencies
io_unit_mod.o:
os_mod.o:               io_unit_mod.o mpi_mod.o omp_mod.o
io_mod.o:               omp_mod.o omp_timer_mod.o mpi_mod.o io_unit_mod.o timer_mod.o os_mod.o
trace_mod.o:            io_mod.o timer_mod.o
timer_mod.o:            io_unit_mod.o omp_mod.o omp_timer_mod.o omp_lock_mod.o io_unit_mod.o mpi_mod.o
legacy_io_mod.o:        io_mod.o trace_mod.o patch_mod.o kinds_mod.o
buffered_io_mod.o:      io_mod.o trace_mod.o patch_mod.o
direct_io_mod.o:        io_mod.o trace_mod.o patch_mod.o
parallel_io_mod.o:      io_mod.o trace_mod.o patch_mod.o mpi_io_mod.o time_slices_mod.o bits_mod.o kinds_mod.o \
                        mpi_file_mod.o mpi_coords_mod.o
data_io_mod.o:          io_mod.o trace_mod.o patch_mod.o os_mod.o omp_mod.o omp_timer_mod.o mpi_mod.o \
                        direct_io_mod.o parallel_io_mod.o buffered_io_mod.o legacy_io_mod.o dll_mod.o \
                        counters_mod.o time_slices_mod.o omp_lock_mod.o pdf_io_mod.o \
                        amr_io_mod.o
validate_mod.o:         io_mod.o io_unit_mod.o patch_mod.o link_mod.o kinds_mod.o
h5a_mod.o:              omp_lock_mod.o
aux_mod.o:              io_unit_mod.o trace_mod.o dll_mod.o
h5_mod.o:               io_mod.o trace_mod.o omp_lock_mod.o
pdf_io_mod.o:           io_mod.o trace_mod.o patch_mod.o counters_mod.o mpi_mod.o \
                        bits_mod.o
used_mod.o:             io_unit_mod.o omp_timer_mod.o
amr_io_mod.o:           io_unit_mod.o trace_mod.o patch_mod.o mpi_io_mod.o \
                        mpi_mod.o list_mod.o time_slices_mod.o
