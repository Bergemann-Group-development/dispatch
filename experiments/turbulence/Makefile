
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
UPD             = adams_bashforth
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
TOP             = ../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o cartesian_mod.o dispatcher_mod.o lagrange_mod.o
experiment_mod.o        : io_mod.o io_unit_mod.o trace_mod.o omp_mod.o bits_mod.o solver_mod.o initial_mod.o \
                          timer_mod.o link_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : forces_mod.o
