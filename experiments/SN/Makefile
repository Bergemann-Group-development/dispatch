
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
UPD             = adams_bashforth
#SOLVER          = stagger2/e
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
TOP             = ../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : dispatcher_mod.o io_mod.o trace_mod.o scaling_mod.o \
                          timer_mod.o cartesian_mod.o SN_mod.o SN_mod.o
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o scaling_mod.o
extras_mod.o            : gpatch_mod.o SN_mod.o cooling_mod.o
SN_mod.o                : gpatch_mod.o scaling_mod.o
sinks_mod.o             : io_mod.o trace_mod.o patch_mod.o scaling_mod.o
scaling_mod.o           : io_mod.o units_mod.o
