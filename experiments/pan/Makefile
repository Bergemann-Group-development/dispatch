
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
#PATCH           = dispatchers
UPD             = adams_bashforth
#SOLVER          = ramses/hydro
SOLVER          = ramses/mhd
#SOLVER          = stagger2/e
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
experiment_mod.o        : io_mod.o trace_mod.o omp_mod.o bits_mod.o solver_mod.o initial_mod.o \
                          timer_mod.o link_mod.o particle_patch_mod.o trace_particles_mod.o validate_mod.o
scaling_mod.o           : io_mod.o units_mod.o
particle_list_mod.o     : omp_timer_mod.o dll_mod.o
trace_particles_mod.o   : io_mod.o trace_mod.o dll_mod.o particle_list_mod.o omp_timer_mod.o random_mod.o mesh_mod.o
