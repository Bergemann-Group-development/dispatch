
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER          = zeus3d
SOLVER          = stagger2
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): components_mod.o setup_mod.o

# Dependencies
components_mod.o:	io_mod.o trace_mod.o rubiks_mod.o
experiment_mod.o	: trace_mod.o io_mod.o kinds_mod.o mesh_mod.o mhd_mod.o pboundary_mod.o
pboundary_mod.o		: io_mod.o trace_mod.o patch_mod.o mesh_mod.o bits_mod.o
initial_mod.o:	io_mod.o trace_mod.o kinds_mod.o mesh_mod.o
