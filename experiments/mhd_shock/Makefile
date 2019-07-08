
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = zeus3d
SOLVER          = stagger2
#SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
#SOLVER          = zeus3d/zeustw
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o	: trace_mod.o io_mod.o kinds_mod.o mesh_mod.o solver_mod.o index_mod.o pboundary_mod.o mpi_mod.o bits_mod.o
pboundary_mod.o		: io_mod.o trace_mod.o mesh_mod.o kinds_mod.o
