
# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
RT_SOLVER       = hybrid_nc
PDE             = solver_mod.o
STAGGER         = stagger_2nd.o
SOLVER          = stagger2
#SOLVER          = stagger2/e
#SOLVER          = AN/hllc
EXEC            = dispatch.x
EOS             = microphysics/eos/ideal
OPACITY         = stub
TOP             = ../..

# include config
sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o eos_mod.o dispatcher_mod.o scaling_mod.o

# Dependencies
initial_mod.o           : io_mod.o trace_mod.o scaling_mod.o mesh_mod.o lagrange_mod.o \
                          atmosphere_mod.o planet_mod.o units_mod.o kinds_mod.o \
                          index_mod.o eos_mod.o
force_mod.o             : io_mod.o trace_mod.o planet_mod.o scaling_mod.o mesh_mod.o
setup_mod.o             : io_mod.o trace_mod.o units_mod.o scaling_mod.o atmosphere_mod.o \
                          rubiks_mod.o planet_mod.o omp_mod.o shared_mod.o list_mod.o \
                          eos_mod.o
atmosphere_mod.o        : io_mod.o trace_mod.o units_mod.o scaling_mod.o particle_mod.o planet_mod.o
experiment_mod.o        : io_mod.o trace_mod.o units_mod.o scaling_mod.o pebble_patch_mod.o kinds_mod.o \
                          timestep_mod.o initial_mod.o atmosphere_mod.o boundaries_mod.o force_mod.o \
                          scalar_mod.o
scaling_mod.o           : io_mod.o units_mod.o math_mod.o
pebbles_mod.o           : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_trace_mod.o      : io_mod.o trace_mod.o bits_mod.o omp_lock_mod.o pebble_mod.o
pebble_patch_mod.o      : io_mod.o trace_mod.o omp_mod.o omp_timer_mod.o task_mod.o $(PDE) link_mod.o planet_mod.o \
                          scaling_mod.o units_mod.o random_mod.o shared_mod.o pebble_mod.o pebbles_mod.o \
                          pebble_trace_mod.o solver_mod.o
pebble_mod.o            :
