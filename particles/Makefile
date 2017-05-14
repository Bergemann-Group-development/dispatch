# Set search path
VPATH := $(VPATH):$(TOP)/particles

# Dependencies
particle_mod.o:         task_mod.o io_mod.o
star_mod.o:             particle_mod.o units_mod.o scaling_mod.o
planet_mod.o:           particle_mod.o units_mod.o scaling_mod.o
pebble_mod.o:           particle_mod.o units_mod.o scaling_mod.o
dust_mod.o:             particle_mod.o
