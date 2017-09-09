
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

RT_PREC=RT_DOUBLE
# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2
#RT_SOLVER       = hybrid
#RT_SOLVER       = hybrid2
RT_SOLVER       = hybrid_nc
#EOS             = microphysics/eos/stagger
#EOS             = microphysics/eos/ideal
EXEC            = dispatch.x
# Ideal or tabular EOS
#EOS             = microphysics/eos/ideal
#OPACITY         = stagger
# Tomida & Hori, with opacity
#EOS              = data/eos/Tomida+Hori_2016/e
EOS             = microphysics/eos/rt_shadow_bench
OPACITY         = stub
#OPACITY         = Tomida+Hori_2016
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
experiment_mod.o	: rmhd_mod.o trace_mod.o io_mod.o initial_mod.o bits_mod.o boundaries_mod.o scaling_mod.o units_mod.o
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o units_mod.o scaling_mod.o
#data_hub_mod.o:		io_mod.o omp_mod.o bits_mod.o kinds_mod.o mesh_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o
eos_mod.o               : io_mod.o mpi_mod.o trace_mod.o scaling_mod.o units_mod.o kinds_mod.o
