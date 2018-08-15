# Append to VPATH
VPATH := $(VPATH):$(TOP)/states

# Dependencies
variable_mod.o:		 mpi_mod.o trace_mod.o kinds_mod.o
state_mod.o:             variable_mod.o io_mod.o trace_mod.o
global_mod.o:		 trace_mod.o lagrange_mod.o
bits_mod:
