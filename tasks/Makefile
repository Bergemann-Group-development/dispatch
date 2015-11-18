# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Local object files to include
ifeq ($(EXEC),dispatch.x)
OBJ += \
  bits_mod.o \
  task_mod.o \
  tasks_mod.o \
  patch_mod.o \
  unigrid_mod.o \
  hd_mod.o \
  mhd_mod.o \
  particle_mod.o \
  dust_mod.o \
  pebble_mod.o \
  planet_mod.o \
  star_mod.o
endif

# Dependencies
item_mod.o:           bits_mod.o trace_mod.o
task_mod.o:           io_mod.o mesh_mod.o bits_mod.o
tasks_mod.o:          task_mod.o
patch_mod.o:          task_mod.o mesh_mod.o
hd_mod.o:             patch_mod.o trace_mod.o
mhd_mod.o:            patch_mod.o trace_mod.o
unigrid_mod.o:        patch_mod.o mpi_coord_mod.o

particle_mod.o:       tasks_mod.o io_mod.o
star_mod.o:           particle_mod.o
planet_mod.o:         particle_mod.o
pebble_mod.o:         particle_mod.o
dust_mod.o:           particle_mod.o
