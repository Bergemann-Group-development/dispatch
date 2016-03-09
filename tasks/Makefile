# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Dependencies
item_mod.o:             bits_mod.o trace_mod.o
task_mod.o:             io_mod.o mesh_mod.o bits_mod.o mpi_mod.o
patch_mod.o:            task_mod.o mesh_mod.o
hd_mod.o:               patch_mod.o trace_mod.o
unigrid_mod.o:          patch_mod.o mpi_coords_mod.o
refine_mod.o:           io_mod.o trace_mod.o patch_mod.o experiment_mod.o
scaling_mod.o:          io_mod.o cgs_mod.o
download_mod.o:         io_mod.o trace_mod.o task_mod.o patch_mod.o list_mod.o
random_mod.o:           io_mod.o trace_mod.o
lagrange_mod.o:         io_mod.o
pebbles_mod.o:          io_mod.o trace_mod.o omp_mod.o omp_timer_mod.o patch_mod.o planet_mod.o star_mod.o \
                        scaling_mod.o cgs_mod.o random_mod.o list_mod.o task_mod.o mhd_mod.o
cgs_mod.o:				
