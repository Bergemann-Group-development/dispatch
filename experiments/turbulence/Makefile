
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
UPD             = adams_bashforth
SOLVER          = ramses/hydro
TOP             = ../..

# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o cartesian_mod.o
experiment_mod.o        : io_mod.o io_unit_mod.o trace_mod.o omp_mod.o bits_mod.o mhd_mod.o initial_mod.o timer_mod.o link_mod.o
