# Append to VPATH, first specific, then general
VPATH := $(VPATH):$(TOP)/solvers/particles/$(PARTICLE_SOLVER):$(TOP)/solvers/particles

# include solver Makefile
sinclude $(TOP)/solvers/particles/$(PARTICLE_SOLVER)/Makefile

# solver info
info::
	@echo "$B  PARTICLE_SOLVER $N= $(PARTICLE_SOLVER)"

# dependencies
patch_force_mod.o       : dll_mod.o
patch_forces_mod.o      : patch_force_mod.o
particle_mod.o          : patch_forces_mod.o io_unit_mod.o lagrange_mod.o
particles_mod.o         : particle_mod.o hash_table_mod.o patch_mod.o
particle_solver_mod.o   : particles_mod.o
