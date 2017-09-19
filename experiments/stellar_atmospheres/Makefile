
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
EXEC            = dispatch.x
TOP             = ../..

# Ideal or tabular EOS
#EOS             = microphysics/eos/ideal
#OPACITY         = stub
# Legacy stagger code
EOS             = microphysics/eos/stagger
OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

# Include Makefile chain
sinclude $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o setup_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : mhd_mod.o trace_mod.o io_mod.o initial_mod.o download_mod.o link_mod.o
initial_mod.o           : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o stagger_mod.o mesh_mod.o
stagger_mod.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o \
                          stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o patch_mod.o scaling_mod.o units_mod.o download_mod.o \
                          stagger_mod.o mesh_mod.o
patches_mod.o           : io_mod.o trace_mod.o timestep_mod.o
eos_mod.o               : io_mod.o trace_mod.o
