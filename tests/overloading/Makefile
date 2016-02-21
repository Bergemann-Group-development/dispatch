# This Makefile is a templater for simple tests, which do not need the complete DISPATCH list of object files

# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized
#OPTS = 

# This should be different from dispatch.x for simpe tests, to disable the normal OBJ list
EXEC = test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

$(EXEC): scalar_mod.o vector_mod.o mesh_mod.o
