# This Makefile is a template for experiments, and may only need to be changed in the Dependencies section

# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized
#OPTS =

# Problem specification
MESH = centered
PDE = advection
TIMESTEP = adams_bashforth

# Common configuration
EXEC = dispatch.x
AUTO = yes
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies
$(EXEC) : experiment_mod.o
experiment_mod.o : io_mod.o task_mod.o patch_mod.o task_list_mod.o trace_mod.o mpi_coords_mod.o
