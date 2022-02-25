
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = stagger2/bifrost
SOLVER          = ramses/hlls
#STAGGER         = stagger_6th.o
RT_SOLVER       = short_characteristics
#RT_SOLVER       = hybrid
EXEC            = dispatch.x
EOS             = microphysics/eos/rt_shadow_bench
OPACITY         = stub
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)            : restart_mod.o cartesian_mod.o setup_mod.o eos_mod.o dispatcher_mod.o

# Other local ependencies
initial_condition_mod.o	: io_mod.o trace_mod.o kinds_mod.o eos_mod.o mesh_mod.o patch_mod.o index_mod.o units_mod.o
extras_mod.o       : patch_mod.o gpatch_mod.o connect_mod.o rt_mod.o initial_condition_mod.o
experiment_mod.o   : link_mod.o solver_mod.o download_mod.o scaling_mod.o eos_mod.o rt_solver_mod.o
scaling_mod.o      : io_mod.o trace_mod.o units_mod.o
boundaries_mod.o          : standard_boundaries_mod.o io_mod.o trace_mod.o \
                            mpi_mod.o bits_mod.o mesh_mod.o index_mod.o \
                            eos_mod.o vector_ops.o units_mod.o scaling_mod.o
