
# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
RT_SOLVER       = short_characteristics
SOLVER          = ramses/hlld_eos
#SOLVER          = stagger2/bifrost
#STAGGER         = stagger_6th.o
EXEC            = dispatch.x
TOP             = ../..
conductivity    = extras/heat_conductivity/spitzer

# Ideal or tabular EOS
#EOS             = microphysics/eos/ideal
#OPACITY         = stub
# Legacy stagger code
EOS             = microphysics/eos/stagger
OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016
CONDUCTIVITY     = spitzer

# Include Makefile chain
sinclude $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile


# List the .o files the main program depends on after the colon:
$(EXEC)              : restart_mod.o cartesian_mod.o setup_mod.o dispatcher_mod.o \
                       read_snapshot_mod.o

# Other local ependencies
extras_mod.o         : patch_mod.o gpatch_mod.o connect_mod.o gravity_mod.o \
                       newton_mod.o rt_mod.o spitzer_mod.o initial_conditions_mod.o \
                       pressure_node_mod.o
boundaries_mod.o     : standard_boundaries_mod.o io_mod.o trace_mod.o \
                       mpi_mod.o bits_mod.o mesh_mod.o index_mod.o \
                       eos_mod.o vector_ops.o units_mod.o scaling_mod.o
experiment_mod.o     : link_mod.o solver_mod.o download_mod.o scaling_mod.o \
                       eos_mod.o rt_solver_mod.o vector_ops.o timestep_mod.o \
                       initial_mod.o
scaling_mod.o        : io_mod.o trace_mod.o units_mod.o
initial_conditions.o : io_mod.o trace_mod.o
