
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = ramses/hydro
EXEC            = dispatch.x
TOP             = ../../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o cartesian_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o timestep_mod.o initial_mod.o download_mod.o link_mod.o

