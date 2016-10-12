
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
SOLVER          = stagger2
RT_SOLVER       = cooling_function
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o timestep_mod.o initial_mod.o download_mod.o link_mod.o 
boundary_mod.o          : patch_mod.o
