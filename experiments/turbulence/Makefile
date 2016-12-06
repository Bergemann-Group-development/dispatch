
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
#EOS             = tabular_e
#OPACITY         = tabular_e
# Tomida & Hori, with opacity
#EOS_DATA        = Tomida+Hori_2016
#OPACITY_DATA    = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : io_mod.o mpi_mod.o cartesian_mod.o
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o initial_mod.o download_mod.o link_mod.o eos_mod.o

