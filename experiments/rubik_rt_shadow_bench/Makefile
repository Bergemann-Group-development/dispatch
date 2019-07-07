
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = stagger2/e
RT_SOLVER       = short_characteristics
#RT_SOLVER       = hybrid
#EOS             = stagger
#EOS             = ideal
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
$(EXEC): setup_mod.o eos_mod.o dispatcher_mod.o rt_setup.o cartesian_mod.o

# Dependencies
extras_mod.o            : patch_mod.o gpatch_mod.o connect_mod.o rt_mod.o
boundary_mod.o          : io_mod.o trace_mod.o units_mod.o scaling_mod.o
setup_mod.o             : io_mod.o trace_mod.o units_mod.o scaling_mod.o omp_mod.o
experiment_mod.o	: link_mod.o solver_mod.o download_mod.o scaling_mod.o
intial_mod.o		: io_mod.o trace_mod.o mesh_mod.o units_mod.o scaling_mod.o

