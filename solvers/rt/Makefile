# Append the rt VPATH, in case there are shared modules
VPATH := $(VPATH):$(TOP)/solvers/rt

# Prepend RT_SOLVER path, allowing sharing of modules
VPATH := $(TOP)/solvers/rt/$(RT_SOLVER):$(VPATH)

# Add to rules for "make info"
RT_PREC ?= RT_SINGLE
info::
	@echo "$B  RT_SOLVER   $N= $(RT_SOLVER)"
	@echo "$B  RT_PREC     $N= $(RT_PREC)"

# Source the RT_SOLVER Makefile
sinclude $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile

# Dependencies on shared RT modules
radau_mod.o         :
