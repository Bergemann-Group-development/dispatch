
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
UPD             = adams_bashforth
TOP             = ../../..

sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o cartesian_mod.o

experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o initial_mod.o pboundary_mod.o
pboundary_mod.o         : io_mod.o trace_mod.o mesh_mod.o patch_mod.o bits_mod.o initial_mod.o
index_mod.o             :
solver_mod.o            : kinds_mod.o gpatch_mod.o
