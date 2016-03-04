# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
stagger_mod.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o patch_mod.o stagger_mod.o boundary_mod.o
patches_mod.o           : io_mod.o trace_mod.o timestep_mod.o
initial_mod.o           : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o
force_mod.o             : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o
boundary_mod.o           : io_mod.o trace_mod.o mesh_mod.o
