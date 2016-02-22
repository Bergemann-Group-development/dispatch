# Append to VPATH
VPATH := $(VPATH):$(TOP)/lists

# Dependencies
list_mod.o        : task_mod.o particle_mod.o planet_mod.o star_mod.o pebble_mod.o patch_mod.o
task_list_mod.o   : mpi_mod.o list_mod.o task_mod.o bits_mod.o trace_mod.o io_mod.o refine_mod.o
mesg_list_mod.o	  : task_list_mod.o
nbor_list_mod.o   : simple_list_mod.o task_list_mod.o task_mod.o
patch_list_mod.o  : experiment_mod.o io_mod.o bits_mod.o patch_mod.o refine_mod.o
