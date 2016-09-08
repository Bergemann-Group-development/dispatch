
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2
MESH            = staggered
EXEC            = dispatch.x

TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o

# Dependencies
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o timestep_mod.o initial_mod.o bits_mod.o
initial_mod.o           : io_mod.o trace_mod.o mesh_mod.o mhd_mod.o cgs_mod.o scaling_mod.o mpi_mod.o
setup_mod.o             : io_mod.o trace_mod.o cgs_mod.o experiment_mod.o task_mod.o task_list_mod.o mesh_mod.o scaling_mod.o bits_mod.o
scaling_mod.o           : io_mod.o trace_mod.o
