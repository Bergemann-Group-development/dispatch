
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized
#OPTS=

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = zeus3d
#SOLVER          = stagger2/bifrost
SOLVER          = ramses/mhd
STAGGER         = stagger_2nd.o
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): restart_mod.o cartesian_mod.o setup_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o	: trace_mod.o io_mod.o kinds_mod.o mesh_mod.o solver_mod.o 
pboundary_mod.o		: io_mod.o trace_mod.o patch_mod.o mesh_mod.o bits_mod.o
initial_mod.o:	io_mod.o trace_mod.o kinds_mod.o mesh_mod.o index_mod.o math_mod.o
