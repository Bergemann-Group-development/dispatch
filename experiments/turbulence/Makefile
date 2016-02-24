
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger
SOLVER          = stagger2
MESH            = staggered
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): task_list_mod.o $(PICOBJ)

# Dependencies
dispatch.o              : task_list_mod.o
experiment_mod.o        : mhd_mod.o scalar_mod.o trace_mod.o vector_mod.o forcing_mod.o io_mod.o timestep_mod.o
forcing_mod.o           : random_mod.o io_mod.o vector_mod.o mesh_mod.o trace_mod.o scalar_mod.o
random_mod.o            : io_mod.o
mhd_mod.o               : patch_mod.o io_mod.o vector_mod.o trace_mod.o stagger_mod.o scalar_mod.o

