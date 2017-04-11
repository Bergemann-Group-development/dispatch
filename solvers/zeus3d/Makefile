# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# **Don't forget to update the LD_LIBRARY_PATH environment variable when you run the code!**
LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore
PREC=DOUBLE

# Dependencies
mhd_mod.o:      io_mod.o trace_mod.o gpatch_mod.o kinds_mod.o variable_mod.o mesh_mod.o timestep_mod.o mpi_mod.o index_mod.o
initial_mod.o:  io_mod.o trace_mod.o variable_mod.o mesh_mod.o random_mod.o index_mod.o
index_mod.o:

