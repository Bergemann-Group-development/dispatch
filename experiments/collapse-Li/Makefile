#COMPILER = gfort10
#COMPILER = ifort

# Option groups
OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
#UPD             = adams_bashforth
#SOLVER          = zeus3d
SOLVER          = ramses/mhd
#SOLVER          = stagger2/e
#SOLVER          = stagger2
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o cartesian_mod.o dispatcher_mod.o initial_conditions_mod.o

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o kinds_mod.o bits_mod.o solver_mod.o pboundary_mod.o
initial_condition_mod.o : io_mod.o trace_mod.o kinds_mod.o mesh_mod.o patch_mod.o index_mod.o units_mod.o
scaling_mod.o           : io_mod.o units_mod.o math_mod.o
pboundary_mod.o         : io_mod.o trace_mod.o mesh_mod.o gpatch_mod.o bits_mod.o
extras_mod.o            : patch_mod.o selfgravity_mod.o scratch_mod.o initial_conditions_mod.o
