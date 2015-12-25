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
$(EXEC): experiment_mod.o
experiment_mod.o : nbor_list_mod.o mhd_mod.o advection_mod.o task_list_mod.o unigrid_mod.o
