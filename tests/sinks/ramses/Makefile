
# Option groups
#OPTS = full_debug
#OPTS = debug
OPTS = optimized

# Default target and methods
EXEC            = test.x
UPD             = adams_bashforth
#SOLVER          = stagger2/e
SOLVER          = ramses/hydro
SINK_PARTICLES  = ramses
#SOLVER          = ramses/mhd
COOLING         = simple
TOP             = ../../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : sinks_mod.o
scaling_mod.o           : io_mod.o units_mod.o
