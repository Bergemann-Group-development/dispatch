# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Local object files to include
#OBJ += \
#task_mod.o \
#hd_mod.o \
#mhd_mod.o \
#list_mod.o \
#patch_mod.o \
#particle_mod.o \
#planet_mod.o \
#star_mod.o \
#mpi_timer.o

# Dependencies
task_mod.f90:           item_mod.o io_mod.o
patch_mod.f90:          task_mod.o mesh_mod.o
hd_mod.f90:             patch_mod.o trace_mod.o
mhd_mod.f90:            patch_mod.o trace_mod.o

particle_mod.f90:       task_mod.o
star_mod.f90:           particle_mod.o
planet_mod.f90:         particle_mod.o
pebble_mod.f90:         particle_mod.o
dust_mod.f90:           particle_mod.o
