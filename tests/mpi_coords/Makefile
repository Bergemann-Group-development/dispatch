# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

default: test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = omp_mod.o mpi_basic.o mpi_coords.o
test.x: $(OBJ)
