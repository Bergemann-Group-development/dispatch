
# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
#SOLVER          = zeus3d
#SOLVER          = stagger2/bifrost
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd

# Standard entries -- do not edit
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o	: solver_mod.o
extras_mod.o            : gpatch_mod.o
