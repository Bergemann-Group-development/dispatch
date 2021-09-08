SINK_PARTICLES ?= ramses

VPATH := $(VPATH):$(TOP)/extras/sink_particles/$(SINK_PARTICLES)

sinclude $(TOP)/extras/sink_particles/$(SINK_PARTICLES)/Makefile

sinks_mod.o: SNe_mod.o patch_force_mod.o particle_mod.o particle_solver_mod.o
