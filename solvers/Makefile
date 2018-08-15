
sinclude $(TOP)/solvers/$(SOLVER)/Makefile
sinclude $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile
sinclude $(TOP)/solvers/poisson/Makefile
sinclude $(TOP)/solvers/particles/Makefile

VPATH := $(VPATH):$(TOP)/solvers

stagger_mod.o           : io_mod.o trace_mod.o stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
stagger_4th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
