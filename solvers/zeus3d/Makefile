# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore
PREC=DOUBLE

# Dependencies
mhd_mod.o:		io_mod.o trace_mod.o patch_mod.o kinds_mod.o variable_mod.o mesh_mod.o download_mod.o
initial_mod.o:	io_mod.o trace_mod.o variable_mod.o mesh_mod.o random_mod.o mhd_mod.o
