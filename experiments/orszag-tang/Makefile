
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger
SOLVER          = stagger2
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o

# Dependencies
experiment_mod.o        : io_mod.o bits_mod.o mhd_mod.o trace_mod.o timestep_mod.o initial_mod.o

