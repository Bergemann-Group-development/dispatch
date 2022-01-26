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

#OPT = -O0 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -fbacktrace

# List the .o files the main program depends on after the colon:
$(EXEC): read_external_mod.o

# Dependencies
read_external_mod.o     : io_mod.o m3io_mod.o mesh_mod.o task_list_mod.o experiment_mod.o \
                          trace_mod.o task_mod.o link_mod.o bits_mod.o halo_mod.o scaling_mod.o \
                          units_mod.o timer_mod.o patch_mod.o download_mod.o list_mod.o \
                          m3d_patch_mod.o
experiment_mod.o        : trace_mod.o io_mod.o solver_mod.o kinds_mod.o atom_var.o spectrum_mod.o \
                          m3d_patch_mod.o opac_patch_mod.o rad_patch_mod.o atom_patch_mod.o rate_patch_mod.o
dispatch.o              : io_mod.o mpi_mod.o
boundaries_mod.o        : standard_boundaries_mod.o io_mod.o trace_mod.o \
                          mpi_mod.o bits_mod.o mesh_mod.o index_mod.o units_mod.o \
                          eos_mod.o vector_ops.o units_mod.o scaling_mod.o
extras_mod.o            : newton_mod.o gravity_mod.o
guards_mod.o            : math_mod.o
coll_var.o              : coll_mod.o
coll_mod.o              : mmath.o 
atom_var.o              : m3io_mod.o spectrum_mod.o trans_mod.o coll_mod.o
spectrum_mod.o          : mmath.o trans_mod.o
m3d_patch_mod.o         : gpatch_mod.o composition_mod.o trace.o
opac_patch_mod.o        : m3d_patch_mod.o
rad_patch_mod.o         : m3d_patch_mod.o geometry_var.o mgeometry.o
rate_patch_mod.o        : m3d_patch_mod.o
atom_patch_mod.o        : rate_patch_mod.o
mmath.o                 : search_mod.o
