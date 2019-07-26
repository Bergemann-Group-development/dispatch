
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
RT_SOLVER       = short_characteristics
#RT_SOLVER       = hybrid
EXEC            = dispatch.x
EOS             = microphysics/eos/rt_shadow_bench
OPACITY         = stub
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)            : cartesian_mod.o setup_mod.o eos_mod.o dispatcher_mod.o rt_nbors_mod.o

# Other local ependencies
initial_mod.o      : io_mod.o trace_mod.o mesh_mod.o units_mod.o scaling_mod.o index_mod.o
extras_mod.o       : patch_mod.o gpatch_mod.o connect_mod.o rt_mod.o
experiment_mod.o   : link_mod.o solver_mod.o download_mod.o scaling_mod.o eos_mod.o rt_solver_mod.o rt_nbors_mod.o
scaling_mod.o      : io_mod.o trace_mod.o units_mod.o
