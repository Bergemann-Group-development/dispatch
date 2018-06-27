
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
EXEC            = dispatch.x
TOP             = ../../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016/e

# Include Makefile chain
sinclude $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): dispatcher_mod.o cartesian_mod.o setup_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : io_mod.o io_unit_mod.o trace_mod.o omp_mod.o bits_mod.o solver_mod.o initial_mod.o \
                          timer_mod.o link_mod.o timestep_mod.o particle_patch_mod.o
scaling_mod.o		: io_unit_mod.o units_mod.o
