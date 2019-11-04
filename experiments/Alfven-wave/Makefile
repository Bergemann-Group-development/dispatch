
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER         = stagger2/e
#SOLVER          = ramses/mhd
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): dispatcher_mod.o setup_mod.o cartesian_mod.o

# Dependencies
experiment_mod.o        : solver_mod.o solver_mod.o trace_mod.o force_mod.o io_mod.o initial_mod.o boundary_mod.o
boundary_mod.o		: io_mod.o trace_mod.o mesh_mod.o solver_mod.o index_mod.o
initial_mod.o           : io_mod.o trace_mod.o force_mod.o mesh_mod.o random_mod.o omp_mod.o index_mod.o
