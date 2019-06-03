
sinclude $(TOP)/solvers/$(SOLVER)/Makefile
sinclude $(TOP)/solvers/particles/Makefile
sinclude $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile
sinclude $(TOP)/solvers/poisson/Makefile
sinclude $(TOP)/solvers/particles/Makefile

VPATH := $(VPATH):$(TOP)/solvers
STAGGER ?= stagger_2nd.o

stagger_2nd.o           : io_mod.o trace_mod.o stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
stagger_4th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
stagger_6th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
vector_ops.o            : $(STAGGER) mesh_mod.o
forces_mod.o            : io_mod.o trace_mod.o patch_mod.o force_mod.o link_mod.o mpi_global_mod.o
