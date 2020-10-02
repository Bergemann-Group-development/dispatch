
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
#PATCH           = dispatchers
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
$(EXEC)                 : setup_mod.o cylinder_mod.o dispatcher_mod.o io_unit_mod.o
experiment_mod.o        : io_mod.o trace_mod.o cpatch_mod.o mesh_mod.o boundaries_mod.o \
                          download_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o list_mod.o forces_mod.o
download_mod.o          : io_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o interpolate_mod.o \
                          kinds_mod.o lagrange_mod.o validate_mod.o
