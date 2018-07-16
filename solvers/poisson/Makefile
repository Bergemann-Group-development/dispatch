# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/poisson

# Dependencies
laplace_mod.o           : io_mod.o trace_mod.o patch_mod.o boundary_mod.o
poisson_mod.o           : io_mod.o trace_mod.o patch_mod.o mesh_mod.o scaling_mod.o
