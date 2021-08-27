
# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and methods
#PATCH           = dispatchers
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
#SOLVER          = stagger2/e
PARTICLE_SOLVER = dust
DUST_EVOLVE     = martin_s
#COOLING         = simple
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
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o mesh_mod.o boundaries_mod.o \
                          vector_mod.o disk_model_mod.o omp_lock_mod.o
disk_model_mod.o        : scaling_mod.o units_mod.o math_mod.o eos_mod.o
disk_particles_mod.o    : disk_model_mod.o vector_ops.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o list_mod.o \
                          selfgravity_mod.o disk_particles_mod.o
#                          forces_mod.o cooling_mod.o
disk_particles_mod.o    : particle_solver_mod.o disk_model_mod.o
disk_selfgravity_mod.o  : selfgravity_mod.o
