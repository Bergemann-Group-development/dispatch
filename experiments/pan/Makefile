
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
#PATCH           = dispatchers
UPD             = adams_bashforth
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
#SOLVER          = stagger2/e
PARTICLE_SOLVER = KDK
EXEC            = dispatch.x
TOP             = ../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o cartesian_mod.o dispatcher_mod.o io_unit_mod.o
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o list_mod.o connect_mod.o forces_mod.o \
                          trace_particles_mod.o sink_patch_mod.o particle_solver_mod.o
