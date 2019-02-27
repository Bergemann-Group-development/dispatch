# Set search path.  Any extras module that needs to have a special 
# version for a solver should be place in solvers/$(SOLVER)/extra,
# and if a local solver dependent version is needed, it should be
# placed in /$(SOLVER)/extra, below the experiment directory
#
VPATH := $(VPATH):$(SOLVER)/extras:$(TOP)/solvers/$(SOLVER)/extras:$(TOP)/extras

# Dependencies
selfgravity_mod.o:      patch_mod.o poisson_mod.o vector_ops.o download_mod.o
