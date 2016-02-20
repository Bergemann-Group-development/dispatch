# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD = adams_bashforth
SOLVERS = compastro
EXEC= dispatch.x

# Common configuration
TOP = ../..
sinclude $(TOP)/config/Makefile

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o patches_mod.o scalar_mod.o ghd_solver_mod.o boundary_mod.o
boundary_mod.o          : trace_mod.o patches_mod.o
dispatch.o              : patch_list_mod.o

# List the .o files the main program depends on after the colon:
$(EXEC): patch_list_mod.o
