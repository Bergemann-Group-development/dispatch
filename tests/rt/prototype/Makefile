
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
EXEC            = dispatch.x
TOP             = ../../..

# Ideal or tabular EOS
#EOS             = ideal
#OPACITY         = stub
# Legacy stagger code
EOS             = stagger
EOS_DATA        = stagger
OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

# Include Makefile chain
sinclude $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o eos_mod.o opacity_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : io_mod.o trace_mod.o timestep_mod.o initial_mod.o rmhd_mod.o
initial_mod.o           : io_mod.o trace_mod.o io_unit_mod.o mesh_mod.o
stagger_mod.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o \
                          stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
rt_solver_mod.o         : io_mod.o trace_mod.o
rt_bundle_mod.o         : io_mod.o trace_mod.o rt_solver_mod.o
rt_mod.o                : io_mod.o trace_mod.o patch_mod.o rt_bundle_mod.o eos_mod.o
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o rt_mod.o

