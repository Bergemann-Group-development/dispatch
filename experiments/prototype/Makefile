# This Makefile is a template for experiments, and may only need to be changed in the Dependencies section

# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized
#OPTS =

# Problem specification
MESH = centered
PROBLEM = advection
TIMESTEP = adams_bashforth

# Local object files
OBJ += experiment_mod.o

# Common configuration
EXEC = dispatch.x
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies
experiment_mod.o : item_mod.o task_mod.o task_list_mod.o trace_mod.o mpi_coord_mod.o
