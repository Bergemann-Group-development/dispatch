# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target
default: objects.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = \
  mpi_timer_mod.o mpi_mod.o mpi_coords_mod.o mpi_send_mod.o mpi_reduce_mod.o \
  upwind_mod.o central_mod.o derivatives_mod.o

# Dependencies
upwind_mod.o: mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o
central_mod.o: mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o

objects.x: $(OBJ)
