
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target and methods
EXEC            = dispatch.x
UPD             = adams_bashforth
SOLVER          = stagger2/e
TOP             = ../..
RT_SOLVER       = hybrid2
# Ideal or tabular EOS
#EOS             = ideal
#OPACITY         = stub
# Tomida & Hori, with opacity
EOS             = tabular_e
EOS_DATA        = Tomida+Hori_2016
OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o microphysics_factory_mod.o scaling_mod.o 
#opacity_mod.o

# Dependencies
#mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o patch_mod.o scaling_mod.o cgs_mod.o download_mod.o \
#                          stagger_mod.o eos_mod.o
microphysics_mod.o      : io_mod.o io_unit_mod.o trace_mod.o scaling_mod.o cgs_mod.o eos_mod.o opacity_mod.o rt_data_mod.o
experiment_mod.o	: mhd_mod.o trace_mod.o io_mod.o initial_mod.o bits_mod.o 
#eos_mod.o opacity_mod.o 
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o mhd_mod.o cgs_mod.o scaling_mod.o


