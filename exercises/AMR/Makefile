# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target
default: experiment.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = \
  mpi_timer_mod.o upwind_mod.o central_mod.o derivatives_mod.o \
  patches_mod.o meshes_mod.o io_mod.o lists_mod.o trace_mod.o omp_mod.o

# Dependencies
lists_mod.o             : io_mod.o patches_mod.o
meshes_mod.o            : trace_mod.o
derivatives_mod.o       : meshes_mod.o
upwind_mod.o            : derivatives_mod.o
central_mod.o           : derivatives_mod.o
patches_mod.o           : derivatives_mod.o upwind_mod.o central_mod.o meshes_mod.o trace_mod.o io_mod.o
experiment.o            : io_mod.o lists_mod.o

experiment.x: $(OBJ)
