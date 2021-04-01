
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#PDE             = mhd_mod.o
STAGGER         = stagger_2nd.o
#STAGGER         = stagger_4th.o
#STAGGER         = stagger_6th.o
SOLVER          = stagger2
#SOLVER          = paper/s
EXEC            = dispatch.x
PARTICLE_SOLVER = KDK
TOP             = ../../..
INITIAL_CONDITIONS = turbulence

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/s
#OPACITY         = Tomida+Hori_2016

# include Makefile chain
include $(TOP)/config/Makefile

#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): dispatcher_mod.o cartesian_mod.o setup_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : io_mod.o io_unit_mod.o trace_mod.o omp_mod.o bits_mod.o solver_mod.o \
                          timer_mod.o link_mod.o omp_timer_mod.o
scaling_mod.o           : units_mod.o
extras_mod.o            : forces_mod.o scratch_mod.o initial_conditions_mod.o
