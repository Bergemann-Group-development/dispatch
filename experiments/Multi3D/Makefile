# Option groups
OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized

# Default target and paths
EXEC            = dispatch.x
TOP             = ../..

# Default methods
SOLVER          = ramses/hlld_eos
RT_SOLVER       = short_characteristics
EOS             = microphysics/eos/stagger
#EXTERNAL        = 1d_model_txt
#
#ifneq ($(EXTERNAL),)
#  VPATH := $(TOP)/components/external/$(EXTERNAL):$(VPATH)
#  sinclude $(TOP)/components/external/$(EXTERNAL)/Makefile
#endif

# To include the Multi3D code, set the macro to point to the src directory
MULTI3D         = $(TOP)/../multi-3d/src
HOP             = ./HOP
sinclude $(TOP)/config/Makefile
sinclude $(HOP)/dispatch.mkf
sinclude $(MULTI3D)/dispatch.mkf

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o read_external_mod.o

# Dependencies
read_external_mod.o     : io_mod.o m3io_mod.o mesh_mod.o task_list_mod.o experiment_mod.o \
						  trace_mod.o task_mod.o link_mod.o bits_mod.o halo_mod.o scaling_mod.o \
						  units_mod.o timer_mod.o patch_mod.o download_mod.o m3d_rt_mod.o list_mod.o
experiment_mod.o        : multi3d_mod.o trace_mod.o io_mod.o solver_mod.o kinds_mod.o
dispatch.o              : io_mod.o mpi_mod.o cartesian_mod.o 
extras_mod.o            : m3d_rt_mod.o newton_mod.o gravity_mod.o
m3d_rt_mod.o            : trace_mod.o rt_mod.o multi3d_mod.o rt_integral_mod.o radau_mod.o math_mod.o \
						  bits_mod.o gpatch_mod.o scs_mod.o list_mod.o link_mod.o halo_mod.o kinds_mod.o
boundaries_mod.o        : standard_boundaries_mod.o io_mod.o trace_mod.o background_mod.o \
                          mpi_mod.o bits_mod.o mesh_mod.o index_mod.o units_mod.o matom.o scs_mod.o \
                          eos_mod.o vector_ops.o units_mod.o scaling_mod.o
guards_mod.o			: math_mod.o
scs_mod.o 				: composition_mod.o
