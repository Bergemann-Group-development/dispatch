
# Option groups
#OPTS = full_debug
OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized
#OPTS = optimized_diagnostics

# TOP directory and default target
TOP             = ../..
EXEC            = dispatch.x

# Choice of solvers
SINK_PARTICLES     = pan
PARTICLE_SOLVER    = KDK
SOLVER             = ramses/hydro
#SOLVER             = ramses/mhd
#SOLVER             = ramses/mhd_eos
#SOLVER             = stagger2
#SOLVER             = stagger2/e
INITIAL_CONDITIONS = turbulence
EOS                = microphysics/eos/ideal

# Tomida & Hori, with opacity
#SOLVER          = ramses/mhd_eos
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

# HDF5 option
#HDF5 = /software/astro/hdf5-1.10.4

sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o restart_mod.o cartesian_mod.o dispatcher_mod.o io_unit_mod.o
experiment_mod.o        : io_mod.o trace_mod.o solver_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o list_mod.o connect_mod.o forces_mod.o \
                          trace_particles_mod.o particle_solver_mod.o SN_mod.o sinks_mod.o selfgravity_mod.o \
                          cooling_mod.o lagrangian_mod.o SN_mod.o initial_conditions_mod.o
particles_mod.o         : particle_mod.o
