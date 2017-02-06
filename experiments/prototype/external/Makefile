
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
EXEC            = dispatch.x
TOP             = ../../..

# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

include $(TOP)/config/Makefile

# Connect to the make system of the code
sinclude code/Makefile
VPATH := $(VPATH):./:./code/

# List the .o files the main program depends on after the colon:
$(EXEC)                 : io_mod.o mpi_mod.o experiment_mod.o cartesian_mod.o
interface_mod.o         :
experiment_mod.o        : io_mod.o trace_mod.o interface_mod.o gpatch_mod.o mpi_mesg_mod.o \
                          dispatch_mod.o io_unit_mod.o timestep_mod.o
