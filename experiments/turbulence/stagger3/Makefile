
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger3/hd
EXEC            = dispatch.x
TOP             = ../../..

# Ideal or tabular EOS
EOS             = ideal
OPACITY         = stub
# Tomida & Hori, with opacity
#EOS             = tabular_e
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016

# include Makefile chain
include $(TOP)/config/Makefile
#PPCODE          = $(TOP)/solvers/pic/PPcode
#sinclude $(PPCODE)/Experiments/DISPATCH/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o eos_mod.o opacity_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : pde_mod.o trace_mod.o io_mod.o timestep_mod.o initial_mod.o download_mod.o link_mod.o scaling_mod.o

