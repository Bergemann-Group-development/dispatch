# default solver
#
PARTICLE_SOLVER ?= void

# Append to specific solver to VPATH
VPATH := $(VPATH):$(TOP)/solvers/particles/$(PARTICLE_SOLVER):$(TOP)/solvers/particles
sinclude          $(TOP)/solvers/particles/$(PARTICLE_SOLVER)/Makefile

# Append generic fall-back
VPATH := $(VPATH):$(TOP)/solvers/particles

# solver info
info::
	@echo "$B  PARTICLE_SOLVER $N= $(PARTICLE_SOLVER)"
