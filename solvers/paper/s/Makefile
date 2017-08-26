# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/paper/Makefile

eos_mod.o       : io_unit_mod.o io_mod.o mpi_mod.o table_mod.o
eos_table_mod.o : io_unit_mod.o io_mod.o mpi_mod.o
rmhd_mod.o      : io_mod.o trace_mod.o mhd_mod.o rt_mod.o
