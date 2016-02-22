
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVERS         = stagger
MESH            = staggered
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
#PPCODE          = $(TOP)/solvers/pic/PPcode
sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): patch_list_mod.o $(PICOBJ)

# Dependencies
dispatch.o              : patch_list_mod.o
experiment_mod.o        : mhd_mod.o scalar_mod.o trace_mod.o vector_mod.o forcing_mod.o io_mod.o
forcing_mod.o           : random_mod.o io_mod.o vector_mod.o mesh_mod.o trace_mod.o scalar_mod.o
random_mod.o            : io_mod.o
refine_mod.o            : trace_mod.o experiment_mod.o io_mod.o patches_mod.o bits_mod.o

# These depend oni -- or are -- local solver/compastro files for now:
mhd_mod.o               : patches_mod.o io_mod.o vector_mod.o trace_mod.o stagger_mod.o scalar_mod.o
#                         -------------
patch_list_mod.o        : io_mod.o patches_mod.o experiment_mod.o trace_mod.o refine_mod.o timestep_mod.o bits_mod.o
#                                  -------------                              ------------
patches_mod.o           : bits_mod.o timestep_mod.o mesh_mod.o omp_timer_mod.o trace_mod.o io_mod.o
# -----------
refine_mod.o            : trace_mod.o experiment_mod.o io_mod.o patches_mod.o bits_mod.o
# ----------                                                    -------------
