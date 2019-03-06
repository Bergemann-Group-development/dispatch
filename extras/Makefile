# Set search path.  Any extras module that needs to have a special 
# version for a solver should be place in solvers/$(SOLVER)/extra,
# and if a local solver dependent version is needed, it should be
# placed in /$(SOLVER)/extra, below the experiment directory
#
VPATH := $(VPATH):$(SOLVER)/extras:$(TOP)/solvers/$(SOLVER)/extras:$(TOP)/extras

# Dependencies
extras_mod.o            : gpatch_mod.o download_mod.o data_hub_mod.o
selfgravity_mod.o       : patch_mod.o poisson_mod.o vector_ops.o download_mod.o
trace_particles_mod.o   : io_mod.o trace_mod.o dll_mod.o particle_list_mod.o omp_timer_mod.o random_mod.o mesh_mod.o link_mod.o patch_mod.o
