
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger
SOLVER          = ramses/mhd
EXEC            = dispatch.x
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): ramses_mod.o download_mod.o

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o scaling_mod.o mhd_mod.o force_mod.o initial_mod.o units_mod.o
ramses_mod.o            : io_mod.o trace_mod.o scaling_mod.o mhd_mod.o force_mod.o initial_mod.o units_mod.o experiment_mod.o star_mod.o task_list_mod.o
scaling_mod.o           : io_mod.o trace_mod.o units_mod.o
