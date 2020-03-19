# Option groups
OPTS = full_debug
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
SOLVER   = ramses/hydro
#SOLVER  = muscl
#SOLVER  = stagger2/e

# Source main configuration Makefile
TOP      = ../../..
sinclude $(TOP)/config/Makefile

# List the .o files the the local files depend on
experiment_mod.o       : io_mod.o trace_mod.o solver_mod.o bits_mod.o mesh_mod.o dust_mod.o
particles_mod.o        : io_mod.o trace_mod.o
particle_solver_mod.o  : io_mod.o trace_mod.o omp_timer_mod.o timer_mod.o task_mod.o \
                         patch_mod.o link_mod.o particles_mod.o
dust_mod.o             : io_mod.o trace_mod.o particle_solver_mod.o
