
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
RT_SOLVER       = hybrid_nc
PDE             = rmhd_mod.o
STAGGER         = stagger_mod.o
SOLVER          = stagger2
#SOLVER          = stagger2/e
#SOLVER          = stagger3/hd/e
#SOLVER          = zeus3d
EXEC            = dispatch.x
EOS             = microphysics/eos/ideal
OPACITY         = stub
TOP             = ../..

# include config
sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o eos_mod.o dispatcher_mod.o scaling_mod.o download_mod.o

# Dependencies
initial_mod.o           : io_mod.o trace_mod.o scaling_mod.o mesh_mod.o lagrange_mod.o atmosphere_mod.o planet_mod.o units_mod.o kinds_mod.o index_mod.o
force_mod.o             : io_mod.o trace_mod.o planet_mod.o scaling_mod.o mesh_mod.o 
setup_mod.o             : io_mod.o trace_mod.o units_mod.o scaling_mod.o atmosphere_mod.o rubiks_mod.o planet_mod.o omp_mod.o scene_mod.o list_mod.o \
                          eos_mod.o radiative_transfer_mod.o
atmosphere_mod.o        : io_mod.o trace_mod.o units_mod.o scaling_mod.o particle_mod.o planet_mod.o
experiment_mod.o        : io_mod.o trace_mod.o units_mod.o scaling_mod.o pebble_patch_mod.o kinds_mod.o \
                          timestep_mod.o initial_mod.o atmosphere_mod.o boundaries_mod.o force_mod.o
scaling_mod.o           : io_mod.o trace_mod.o units_mod.o
units_mod.o               : io_mod.o
star_mod.o              : io_mod.o trace_mod.o units_mod.o scaling_mod.o
planet_mod.o            : io_mod.o trace_mod.o units_mod.o scaling_mod.o
gpatch_mod.o:           io_mod.o trace_mod.o patch_mod.o download_mod.o initial_mod.o kinds_mod.o data_hub_mod.o
pebbles_mod.o           : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_trace_mod.o      : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_patch_mod.o      : io_mod.o trace_mod.o omp_mod.o omp_timer_mod.o task_mod.o $(PDE) link_mod.o planet_mod.o \
                          star_mod.o scaling_mod.o units_mod.o random_mod.o scene_mod.o pebble_mod.o pebbles_mod.o \
                          pebble_trace_mod.o solver_mod.o
pebble_mod.o            :
rt_boundaries_mod.o     : io_mod.o trace_mod.o mpi_mod.o bits_mod.o mesh_mod.o ray_geometry_v3_mod.o planet_mod.o \
                          atmosphere_mod.o eos_mod.o rt_components_mod.o
