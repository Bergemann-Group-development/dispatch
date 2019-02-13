# This Makefile is a templater for simple tests, which do not need the complete DISPATCH list of object files

# Option groups
#OOPTS = debug
#OPTS = full_debug
OPTS = optimized
#OPTS = 

# This should be different from dispatch.x for simpe tests, to disable the normal OBJ list
STAGGER = stagger_2nd.o
EXEC = test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

vector_mod.o: scalar_mod.o trace_mod.o
stagger_mod.o: vector_mod.o trace_mod.o
pde_mod.o: $(STAGGER) trace_mod.o

$(EXEC): mesh_mod.o $(STAGGER) pde_mod.o trace_mod.o
