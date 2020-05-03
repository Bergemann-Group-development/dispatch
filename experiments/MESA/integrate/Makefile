
# Option groups
#OPTS = full_debug
OPTS = debug
#OPTS =
#OPTS = optimized
#OPTS = optimized_diagnostics

# Default target and methods
UPD             = adams_bashforth
#PATCH           = dispatchers
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
#SOLVER          = stagger2/e
PARTICLE_SOLVER = KDK
EXEC            = dispatch.x
TOP             = ../../..

#HDF5 = /software/astro/hdf5-1.10.4
HDF5 =

sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o restart_mod.o cartesian_mod.o dispatcher_mod.o io_unit_mod.o
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o list_mod.o forces_mod.o \
                          selfgravity_mod.o cooling_mod.o mesa_mod.o
sinks_mod.o             : io_mod.o trace_mod.o gpatch_mod.o units_mod.o scaling_mod.o math_mod.o
mesa_mod.o              : sinks_mod.o