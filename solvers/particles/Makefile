# Append to VPATH, first specific, then general
VPATH := $(VPATH):$(TOP)/solvers/particles/$(PARTICLE_SOLVER):$(TOP)/solvers/particles

# include solver Makefile
sinclude $(TOP)/solvers/particles/$(PARTICLE_SOLVER)/Makefile

# solver info
info::
	@echo "$B  PARTICLE_SOLVER $N= $(PARTICLE_SOLVER)"
