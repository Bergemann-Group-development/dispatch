# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/particles

# Dependencies
particle_history_mod.o	: io_mod.o trace_mod.o particle_mod.o patch_mod.o
particle_patch_mod.o	: io_mod.o trace_mod.o particle_history_mod.o particle_mod.o
sinkparticle_mod.o	: io_mod.o trace_mod.o particle_mod.o particle_patch_mod.o patch_mod.o \
                          download_mod.o data_hub_mod.o
