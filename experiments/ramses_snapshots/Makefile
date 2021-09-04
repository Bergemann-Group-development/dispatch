
# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger2
SOLVER          = ramses/mhd
EXTERNAL        = ramses
EXEC            = dispatch.x
TOP             = ../..

sinclude $(TOP)/config/Makefile
VPATH := $(TOP)/particles/:$(VPATH)

# List the .o files the main program depends on after the colon:
$(EXEC): setup_mod.o ramses_mod.o dispatcher_mod.o

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o shared_mod.o bits_mod.o
scaling_mod.o           : io_mod.o trace_mod.o units_mod.o
restart_mod.o           : io_mod.o trace_mod.o task_list_mod.o task_factory_mod.o \
                          amr_io_mod.o read_external_mod.o
