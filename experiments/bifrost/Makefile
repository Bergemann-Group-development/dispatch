
# Option groups
#OPTS = quick
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
SOLVER          = bifrost
UPD             = adams_bashforth
EXEC            = dispatch.x
TOP             = ../..
MESH            = bifrost_mesh

# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

include $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : dispatch_mod.o extern_mod.o cartesian_mod.o $(OBJ) dispatcher_mod.o
dispatch_mod.o		: io_mod.o mpi_mod.o
experiment_mod.o        : trace_mod.o timestep_mod.o solver_mod.o
refine_mod.o:           io_mod.o trace_mod.o bits_mod.o list_mod.o link_mod.o kinds_mod.o \
                        download_mod.o mesh_mod.o solver_mod.o patch_mod.o
