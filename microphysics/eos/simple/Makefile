DATA = $(TOP)/data/eos/simple

sinclude $(DATA)/Makefile
VPATH := $(VPATH):$(DATA)

eos_mod.o       : mpi_mod.o io_unit_mod.o trace_mod.o eos_reader_mod.o scaling_mod.o units_mod.o lagrange_mod.o 
