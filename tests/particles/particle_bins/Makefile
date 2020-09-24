# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

OBJ += list_mod.o task_mod.o particle_mod.o star_mod.o planet_mod.o mpi_timer.o
default: test.x

# Common configuration
TOP = ../../..
include $(TOP)/config/Makefile

test.x: $(OBJ)
