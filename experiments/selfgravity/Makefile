
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
MESH            = staggered
SOLVERS         = compastro
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o patch_mod.o scalar_mod.o ghd_solver_mod.o boundary_mod.o
boundary_mod.o          : trace_mod.o patch_mod.o
dispatch.o              : task_list_mod.o

# List the .o files the main program depends on after the colon:
$(EXEC): task_list_mod.o $(PICOBJ)
