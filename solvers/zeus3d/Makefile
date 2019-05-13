# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/zeus3d

PREC=DOUBLE

# Dependencies
mhd_mod.o:      io_mod.o trace_mod.o force_mod.o extras_mod.o kinds_mod.o variable_mod.o mesh_mod.o timestep_mod.o mpi_mod.o index_mod.o
initial_mod.o:  io_mod.o trace_mod.o variable_mod.o mesh_mod.o random_mod.o index_mod.o
force_mod.o:    io_mod.o trace_mod.o mesh_mod.o random_mod.o
solver_mod.o:   io_mod.o trace_mod.o patch_mod.o mhd_mod.o kinds_mod.o mesh_mod.o task_mod.o
vector_ops.o:	trace_mod.o kinds_mod.o mesh_mod.o

# **Don't forget to update the LD_LIBRARY_PATH environment variable when you run the code!**
ifeq ($(SOLVER),zeus3d)
  LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore
  mhd_mod.o:      FFLAGS+=-I/home/jramsey/repos/ZEUScore
endif
