# This Makefile is a template for experiments, and may only need to be changed in the Dependencies section

# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized
OPTS =

# Problem specification
MESH = centered
PROBLEM = advection
UPD = adams_bashforth
EXEC = dispatch.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies
experiment_mod.o : mpi_mod.o mpi_coord_mod.o list_mod.o mhd_mod.o advection_mod.o task_list_mod.o

$(EXEC): experiment_mod.o
