# Set the RT_SOLVER
#
RT_SOLVER ?= void

VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)
sinclude          $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile

# Append RT_SOLVER path, allowing sharing of modules
VPATH := $(VPATH)::$(TOP)/solvers/rt

# Add to rules for "make info"
RT_PREC ?= RT_SINGLE
info::
	@echo "$B  RT_SOLVER   $N= $(RT_SOLVER)"
	@echo "$B  RT_PREC     $N= $(RT_PREC)"

# Dependencies on shared RT modules
radau_mod.o         :
