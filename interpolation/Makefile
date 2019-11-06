# Set search path
VPATH := $(VPATH):$(TOP)/interpolation

# Dependencies
lagrange_mod.o:         io_mod.o trace_mod.o
interpolate_mod.o:      io_mod.o trace_mod.o mpi_mod.o kinds_mod.o index_mod.o
download_mod.o:         io_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o \
                        interpolate_mod.o kinds_mod.o lagrange_mod.o validate_mod.o \
                        shared_mod.o guard_zones_mod.o remesh_mod.o
data_hub_mod.o:         io_mod.o omp_mod.o bits_mod.o kinds_mod.o mesh_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o
central_mod.o:          patch_mod.o lagrange_mod.o kinds_mod.o
time_slices_mod.o:      io_mod.o omp_mod.o trace_mod.o patch_mod.o kinds_mod.o lagrange_mod.o
remesh_mod.o:           io_mod.o trace_mod.o omp_timer_mod.o patch_mod.o
guard_zones_mod.o:      io_mod.o trace_mod.o omp_timer_mod.o patch_mod.o omp_timer_mod.o
