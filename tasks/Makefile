# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Dependencies
item_mod.o:             bits_mod.o trace_mod.o
task_mod.o:             io_mod.o mesh_mod.o bits_mod.o mpi_mod.o
tasks_mod.o:            task_mod.o
patch_mod.o:            tasks_mod.o mesh_mod.o
hd_mod.o:               patch_mod.o trace_mod.o
mhd_mod.o:              patch_mod.o trace_mod.o
unigrid_mod.o:          patch_mod.o mpi_coords_mod.o
refine_mod.o:           io_mod.o trace_mod.o experiment_mod.o
particle_mod.o:         tasks_mod.o io_mod.o
star_mod.o:             particle_mod.o
planet_mod.o:           particle_mod.o
pebble_mod.o:           particle_mod.o
dust_mod.o:             particle_mod.o
scaling_mod.o:          io_mod.o cgs_mod.o
