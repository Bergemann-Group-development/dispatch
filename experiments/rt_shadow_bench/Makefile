
# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2
RT_SOLVER       = hybrid2
#RT_SOLVER       = hybrid
MESH            = staggered
EXEC            = dispatch.x
# opacities: 
# grey: bell_lin_1994, semenov_2003, rt_shadow_bench
# multi-frequency: vaytet_2015
OPACITIES       = rt_shadow_bench
TOP             = ../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o radiative_transfer_mod.o

# Dependencies
experiment_mod.o	: mhd_rt_mod.o trace_mod.o io_mod.o timestep_mod.o initial_mod.o bits_mod.o
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o mhd_mod.o cgs_mod.o scaling_mod.o

