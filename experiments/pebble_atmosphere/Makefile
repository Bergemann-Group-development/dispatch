
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
SOLVER          = stagger2/e
#SOLVER          = stagger3/hd
EXEC            = dispatch.x
EOS             = ideal
OPACITY         = stub
TOP             = ../..

# include config
sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o eos_mod.o

# Dependencies
setup_mod.o             : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o atmosphere_mod.o rubiks_mod.o planet_mod.o omp_mod.o scene_mod.o list_mod.o
atmosphere_mod.o        : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o particle_mod.o planet_mod.o
experiment_mod.o        : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o pebble_patch_mod.o \
                          timestep_mod.o initial_mod.o atmosphere_mod.o boundary_mod.o force_mod.o
scaling_mod.o           : io_mod.o trace_mod.o cgs_mod.o
cgs_mod.o               : io_mod.o
star_mod.o              : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o
planet_mod.o            : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o
boundary_mod.o          : io_mod.o trace_mod.o cgs_mod.o scaling_mod.o atmosphere_mod.o
pebbles_mod.o           : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_trace_mod.o      : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_patch_mod.o      : io_mod.o trace_mod.o omp_mod.o omp_timer_mod.o task_mod.o mhd_mod.o link_mod.o planet_mod.o \
                          star_mod.o scaling_mod.o cgs_mod.o random_mod.o scene_mod.o pebble_mod.o pebbles_mod.o pebble_trace_mod.o
pebble_mod.o            : 
