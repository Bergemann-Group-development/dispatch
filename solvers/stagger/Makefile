# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
stagger_mod.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o

patches_mod.o           : io_mod.o trace_mod.o timestep_mod.o
