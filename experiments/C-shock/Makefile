
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger
SOLVER          = stagger2
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : mhd_mod.o scalar_mod.o trace_mod.o vector_mod.o force_mod.o io_mod.o timestep_mod.o initial_mod.o boundary_mod.o

