
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
SOLVER = zeus3d
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o cartesian_mod.o

# Dependencies
experiment_mod.o	: io_mod.o trace_mod.o kinds_mod.o mesh_mod.o solver_mod.o
initial_mod.o		: io_mod.o trace_mod.o kinds_mod.o mesh_mod.o patch_mod.o index_mod.o units_mod.o
scaling_mod.o           : io_mod.o units_mod.o
