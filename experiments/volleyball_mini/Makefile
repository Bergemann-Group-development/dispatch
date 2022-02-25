
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized
#OPTS = optimized_diagnostics
OPTS=

# Default target and methods
#PATCH           = dispatchers
UPD             = adams_bashforth
#SOLVER          = ramses/hlld_eos
SOLVER          = ramses/hlld_eos
#SOLVER          = stagger2/e
UPD             = adams_bashforth
#PARTICLE_SOLVER = KDK
EXEC            = dispatch.x
TOP             = ../..
#HDF5 = /software/astro/hdf5-1.10.4
HDF5 =

# Tomida & Hori, with opacity
EOS             = microphysics/eos/itable
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile
VPATH := $(SOLVER):.:$(VPATH)
sinclude $(SOLVER)/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)           : setup_mod.o restart_mod.o cartesian_mod.o dispatcher_mod.o \
                    io_unit_mod.o
cartesian_mod.o         : io_mod.o trace_mod.o experiment_mod.o  task_list_mod.o \
                          mpi_coords_mod.o omp_mod.o scaling_mod.o units_mod.o
scaling_mod.o     : io_mod.o units_mod.o
model_1d_mod.o    : io_mod.o trace_mod.o scaling_mod.o
forces_mod.o      : io_mod.o model_1d_mod.o scratch_mod.o gpatch_mod.o interpolate_mod.o
extras_mod.o      : forces_mod.o newton_mod.o forces_mod.o sphere_point_mod.o \
                    gravity_mod.o rad_diffusion_mod.o
experiment_mod.o  : io_mod.o trace_mod.o solver_mod.o link_mod.o \
                    kinds_mod.o region_mod.o interpolate_mod.o
boundaries_mod.o     : standard_boundaries_mod.o io_mod.o trace_mod.o \
                       mpi_mod.o bits_mod.o mesh_mod.o index_mod.o \
                       eos_mod.o units_mod.o scaling_mod.o model_1d_mod.o interpolate_mod.o
