# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and paths
EXEC            = dispatch.x
TOP             = ../..

# Default methods
SOLVER          = ramses/hlld_eos
RT_SOLVER       = short_characteristics
EOS             = microphysics/eos/stagger

# To include the Multi3D code, set the macro to point to the src directory
MULTI3D         = $(TOP)/../multi-3d/src
sinclude $(TOP)/config/Makefile
sinclude $(MULTI3D)/dispatch.mkf

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o

# Dependencies
experiment_mod.o        : multi3d_mod.o trace_mod.o io_mod.o solver_mod.o
dispatch.o              : io_mod.o mpi_mod.o cartesian_mod.o
extras_mod.o            : m3d_rt_mod.o newton_mod.o gravity_mod.o
m3d_rt_mod.o            : trace_mod.o rt_mod.o multi3d_mod.o
boundaries_mod.o        : vector_ops.o
