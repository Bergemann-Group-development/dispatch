
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized
#OPTS=

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger2/bifrost
SOLVER          = ramses/mhd_eos
STAGGER         = stagger_2nd.o
#STAGGER         = stagger_6th.o
#STAGGER         = stagger_mod.o
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): restart_mod.o cartesian_mod.o setup_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o        : io_mod.o bits_mod.o solver_mod.o trace_mod.o
initial_condition_mod.o	: io_mod.o trace_mod.o kinds_mod.o mesh_mod.o patch_mod.o index_mod.o units_mod.o
extras_mod.o            : patch_mod.o scratch_mod.o initial_conditions_mod.o

