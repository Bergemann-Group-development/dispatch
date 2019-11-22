
# Option groups
OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
RT_SOLVER       = short_characteristics
SOLVER          = stagger2/e
STAGGER         = stagger_2nd.o
EXEC            = dispatch.x
TOP             = ../..

# Ideal or tabular EOS
#EOS             = microphysics/eos/ideal
#OPACITY         = stub
# Legacy stagger code
EOS             = microphysics/eos/stagger
OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

# Include Makefile chain
sinclude $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile


# List the .o files the main program depends on after the colon:
$(EXEC)            : cartesian_mod.o setup_mod.o eos_mod.o dispatcher_mod.o \
                     read_snapshot_mod.o rt_nbors_mod.o

# Other local ependencies
extras_mod.o       : patch_mod.o gpatch_mod.o connect_mod.o gravity_mod.o \
                     newton_mod.o rt_mod.o
experiment_mod.o   : link_mod.o solver_mod.o download_mod.o scaling_mod.o \
                     eos_mod.o rt_solver_mod.o rt_nbors_mod.o
scaling_mod.o      : io_mod.o trace_mod.o units_mod.o
