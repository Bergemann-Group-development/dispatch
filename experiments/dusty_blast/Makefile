
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/multifluid
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o        : solver_mod.o trace_mod.o force_mod.o io_mod.o initial_mod.o timestep_mod.o pboundary_mod.o
force_mod.o 		: io_mod.o trace_mod.o vector_mod.o mesh_mod.o index_mod.o vector_ops.o scalar_mod.o
pboundary_mod.o		: patch_mod.o
