
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = zeus3d
#SOLVER          = stagger2
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o

# Dependencies
experiment_mod.o	: trace_mod.o io_mod.o timestep_mod.o mhd_mod.o initial_mod.o boundary_mod.o
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o mhd_mod.o kinds_mod.o
boundary_mod.o		: io_mod.o trace_mod.o patch_mod.o
