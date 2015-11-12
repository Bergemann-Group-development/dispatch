# This Makefile is a template for experiments, and may only need to be changed in the Dependencies section

# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized
OPTS = 

# Local object files
OBJ += experiment_mod.o
MHD = stagger

# Common configuration
EXEC = dispatch.x
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies
experiment_mod.o : list_mod.o mhd_mod.o
