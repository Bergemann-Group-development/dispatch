# VERBOSE info
ifneq ($(VERBOSE),)
  $(info $(TOP)/solvers/Makefile OFFLOAD=$(OFFLOAD))
endif

# First set VPATH, then allow the sinclude to modify
#
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)
sinclude          $(TOP)/solvers/$(SOLVER)/Makefile

VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)/$(OFFLOAD)
sinclude          $(TOP)/solvers/$(SOLVER)/$(OFFLOAD)/Makefile

# These sub-solvers to the same, internally
sinclude $(TOP)/solvers/rt/Makefile
sinclude $(TOP)/solvers/particles/Makefile
sinclude $(TOP)/solvers/poisson/Makefile

# Allow the top directory to be fall-back
VPATH := $(VPATH):$(TOP)/solvers

STAGGER ?= stagger_2nd.o

stagger_2nd.o           : io_mod.o trace_mod.o stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
stagger_4th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
stagger_6th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
vector_ops.o            : $(STAGGER) mesh_mod.o
forces_mod.o            : io_mod.o trace_mod.o gpatch_mod.o force_mod.o link_mod.o \
                          mpi_global_mod.o scratch_mod.o
flux_mod.o              : io_mod.o trace_mod.o patch_mod.o
divb_clean_mod.o        : io_mod.o trace_mod.o mesh_mod.o vector_mod.o scalar_mod.o vector_ops.o

solver_mod.o            : trace_mod.o # extras_mod.o
divb_cleaner_mod.o      : patch_mod.o stagger_2nd.o kinds_mod.o poisson_mod.o \
                          list_mod.o timer_mod.o