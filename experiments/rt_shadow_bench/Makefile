
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
#RT_SOLVER       = hybrid
#RT_SOLVER       = hybrid2
RT_SOLVER       = hybrid
#EOS             = stagger
#EOS             = ideal
EXEC            = dispatch.x
# Ideal or tabular EOS
#EOS             = ideal
#OPACITY         = stagger
# Tomida & Hori, with opacity
#EOS              = tabular_e
EOS             = rt_shadow_bench
OPACITY         = stub
#OPACITY         = Tomida+Hori_2016
#EOS_DATA        = Tomida+Hori_2016
# opacities: 
# grey: bell_lin_1994, semenov_2003, rt_shadow_bench
# multi-frequency: vaytet_2015
#OPACITIES       = rt_shadow_bench
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o radiative_transfer_mod.o eos_mod.o opacity_mod.o

# Dependencies
boundary_mod.o          : io_mod.o trace_mod.o units_mod.o scaling_mod.o
mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o patch_mod.o scaling_mod.o units_mod.o download_mod.o \
                          stagger_mod.o eos_mod.o
experiment_mod.o	: rmhd_mod.o trace_mod.o io_mod.o initial_mod.o bits_mod.o boundaries_mod.o
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o units_mod.o scaling_mod.o



