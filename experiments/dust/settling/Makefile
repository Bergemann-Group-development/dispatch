# Option groups
#OPTS = full_debug
OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
SOLVER   = ramses/hydro
#SOLVER  = muscl
#SOLVER  = stagger2/e
PARTICLE_SOLVER = dust

# Source main configuration Makefile
TOP      = ../../..
sinclude $(TOP)/config/Makefile

# List the .o files the the local files depend on
experiment_mod.o       : io_mod.o trace_mod.o solver_mod.o bits_mod.o \
                         mesh_mod.o extras_mod.o omp_lock_mod.o
particles_mod.o        : io_mod.o trace_mod.o
extras_mod.o           : gpatch_mod.o particles_mod.o particle_solver_mod.o \
                         list_mod.o
scaling_mod.o          : io_mod.o trace_mod.o math_mod.o units_mod.o
