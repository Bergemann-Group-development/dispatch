
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
SOLVER          = ramses/hydro
TOP             = ../..
EOS_TYPE        = ideal
#EOS_TYPE        = E
#EOS_DATA        = simple
#OPACITY         = semenov_etal_2003
EXEC            = dispatch.x
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : io_mod.o mpi_mod.o cartesian_mod.o
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o initial_mod.o download_mod.o link_mod.o eos_mod.o

