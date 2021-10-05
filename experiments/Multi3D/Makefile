# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
SOLVER          = ramses/hlld_eos
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# To include the Multi3D code, set the macro to point to the src directory
MULTI3D         = $(TOP)/../multi-3d/src
sinclude $(MULTI3D)/dispatch.mkf

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o

# Dependencies
experiment_mod.o        : multi3d_mod.o trace_mod.o io_mod.o solver_mod.o
dispatch.o              : io_mod.o mpi_mod.o cartesian_mod.o

