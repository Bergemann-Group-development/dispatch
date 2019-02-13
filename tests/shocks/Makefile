
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2
#SOLVER          = stagger2/e
#SOLVER          = ramses/mhd
STAGGER         = stagger_2nd.o
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): dispatcher_mod.o setup_mod.o cartesian_mod.o

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o mesh_mod.o bits_mod.o $(STAGGER)
