
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized
#OPTS = optimized_diagnostics

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/bifrost
STAGGER         = stagger_6th.o
MESH            = staggered
EXEC            = dispatch.x
EOS             = microphysics/eos/ideal
OPACITY         = stub
TOP             = ../../..

# sinclude $(TOP)/config/Makefile,
# which sincludes $(TOP)/solvers/Makefile,
# which sincludes $(TOP)/solvers/$(SOLVER)/Makefile,
sinclude $(TOP)/config/Makefile

# Choose PIC components (need to be defined before sinclude)
STAGGER_ORDER   = 6th
INTERPOLATION   = tsc
PPMPI           = mpi
DIM             = 3d
BOUNDARY        = mhd_boundary

# pass to the PPcode make system
sinclude $(TOP)/solvers/pic/Makefile

# List the .o files the main program depends on after the colon: $(PICOBJ)
$(EXEC): cartesian_mod.o setup_mod.o cartesian_mod.o dispatcher_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : trace_mod.o io_mod.o timestep_mod.o solver_mod.o pboundary_mod.o
pboundary_mod.o         : io_mod.o trace_mod.o patch_mod.o mesh_mod.o bits_mod.o
initial_mod.o           : io_mod.o trace_mod.o mesh_mod.o math_mod.o index_mod.o scaling_mod.o vector_ops.o
extras_mod.o            : patch_mod.o gpatch_mod.o connect_mod.o pic_mod.o \
                          pic_particles_mod.o scratch_mod.o
