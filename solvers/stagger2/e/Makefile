# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/stagger2/Makefile

eos_mod.o: io_unit_mod.o io_mod.o mpi_mod.o
