# This Makefile is a templater for simple tests, which do not need the complete DISPATCH list of object files

# Option groups
#OOPTS = debug
#OPTS = full_debug
OPTS = optimized
#OPTS = 

# This should be different from dispatch.x for simpe tests, to disable the normal OBJ list
MESH = staggered
EXEC = test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

$(EXEC): mpi_mod.o io_mod.o
