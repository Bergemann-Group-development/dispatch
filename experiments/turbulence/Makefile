
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
SOLVER          = ramses/hydro
TOP             = ../..
# Ideal or tabular EOS
#EOS             = ideal
#OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : io_mod.o mpi_mod.o cartesian_mod.o scaling_mod.o opacity_mod.o
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o initial_mod.o download_mod.o link_mod.o eos_mod.o opacity_mod.o

