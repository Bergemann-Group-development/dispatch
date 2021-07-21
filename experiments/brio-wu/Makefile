
# Option groups
#OPTS = debug
OPTS =
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = zeus3d
#SOLVER          = ramses/hydro
#SOLVER          = stagger2/bifrost
STAGGER         = stagger_2nd.o
SOLVER          = ramses/mhd_eos/s
#SOLVER          = ramses/mhd
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o dispatcher_mod.o initial_conditions_mod.o

# Dependencies
experiment_mod.o	: trace_mod.o io_mod.o solver_mod.o bits_mod.o mesh_mod.o
initial_mod.o           : io_mod.o trace_mod.o mesh_mod.o index_mod.o math_mod.o kinds_mod.o eos_mod.o
