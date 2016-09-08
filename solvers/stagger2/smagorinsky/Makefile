# include "parent" Makefil
include $(TOP)/solvers/stagger2/e/Makefile

# Append to VPATH
VPATH := $(TOP)/solvers/smagorinsky:$(VPATH)

eos_mod.o: io_unit_mod.o io_mod.o mpi_mod.o
