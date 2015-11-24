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
  mpi_timer_mod.o mpi_mod.o mpi_coords_mod.o mpi_send_mod.o mpi_reduce_mod.o \
  upwind_mod.o central_mod.o derivatives_mod.o patch_mod.o mesh_mod.o io_mod.o trace_mod.o omp_mod.o

# Dependencies
derivatives_mod.o       : mesh_mod.o
upwind_mod.o            : mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o
central_mod.o           : mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o
patch_mod.o             : derivatives_mod.o upwind_mod.o central_mod.o mesh_mod.o

experiment.x: $(OBJ)
