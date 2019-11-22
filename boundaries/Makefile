# Set search path
VPATH := $(VPATH):$(TOP)/boundaries

# Dependencies
boundaries_mod.o        : io_mod.o trace_mod.o mpi_mod.o bits_mod.o mesh_mod.o index_mod.o boundary_wavekill_mod.o 
boundary_wavekill_mod.o : kinds_mod.o
