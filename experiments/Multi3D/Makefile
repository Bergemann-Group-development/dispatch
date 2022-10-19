# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
OPTS =
#OPTS = hoppe
#OPTS = optimized

# Default target and paths
EXEC            = dispatch.x
TOP             = ../..

# Default methods
SOLVER          = ramses/hlld_eos
RT_SOLVER       = short_characteristics
#EOS             = microphysics/eos/itable
EOS             = microphysics/eos/stagger
#EXTERNAL        = 1d_model_txt
#
#ifneq ($(EXTERNAL),)
#  VPATH := $(TOP)/components/external/$(EXTERNAL):$(VPATH)
#  sinclude $(TOP)/components/external/$(EXTERNAL)/Makefile
#endif

# To include the Multi3D code, set the macro to point to the src directory
HOP             = ./HOP
sinclude $(TOP)/config/Makefile
sinclude $(HOP)/dispatch.mkf

#OPT = -O0 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -fbacktrace

# List the .o files the main program depends on after the colon:
$(EXEC): read_external_mod.o mpi_debug_mod.o

# Dependencies
read_external_mod.o         : io_mod.o m3io_mod.o mesh_mod.o task_list_mod.o experiment_mod.o \
                              trace_mod.o task_mod.o link_mod.o bits_mod.o halo_mod.o scaling_mod.o \
                              units_mod.o timer_mod.o patch_mod.o download_mod.o list_mod.o \
                              ray_patch_mod.o read_stagger_snapshot_mod.o
read_stagger_snapshot_mod.o : io_mod.o os_mod.o eos_mod.o scaling_mod.o
experiment_mod.o            : trace_mod.o io_mod.o solver_mod.o kinds_mod.o m3d_halo_mod.o
dispatch.o                  : io_mod.o
boundaries_mod.o            : standard_boundaries_mod.o io_mod.o trace_mod.o \
                              mpi_mod.o bits_mod.o mesh_mod.o index_mod.o units_mod.o \
                              eos_mod.o vector_ops.o units_mod.o scaling_mod.o
extras_mod.o                : newton_mod.o gravity_mod.o composition_patch_mod.o
guards_mod.o                : math_mod.o
coll_var.o                  : coll_mod.o
coll_mod.o                  : mmath.o 
trans_mod.o                 : bbtrans_mod.o coll_mod.o
atom_mod.o                  : m3io_mod.o trans_mod.o coll_mod.o op_math_mod.o spectrum_mod.o
spectrum_mod.o              : io_mod.o timer_mod.o trace_mod.o mconst.o
m3d_halo_mod.o              : halo_mod.o guards_mod.o gpatch_mod.o
m3d_patch_mod.o             : gpatch_mod.o atom_mod.o ray_mod.o m3d_halo_mod.o
ray_patch_mod.o             : m3d_patch_mod.o
rad_patch_mod.o             : ray_patch_mod.o back_opac_patch_mod.o
rate_patch_mod.o            : ray_patch_mod.o rad_patch_mod.o lines_patch_mod.o mean_rad_patch_mod.o
atom_patch_mod.o            : rate_patch_mod.o lines_patch_mod.o cont_opac_patch_mod.o
composition_patch_mod.o     : out_rad_patch_mod.o back_opac_patch_mod.o rad_patch_mod.o atom_patch_mod.o net_opac_patch_mod.o proj_patch_mod.o mean_rad_patch_mod.o
bftrans_mod.o               : os_mod.o
bbtrans_mod.o               : trace_mod.o
mmath.o                     : search_mod.o
unista_mod.o                : refsor_mod.o
ray_mod.o                   : mmath.o
prec_mod.o                  : mconst.o os_mod.o
patch_mod.o                 : os_mod.o
mpi_debug_mod.o             : mpi_mod.o io_mod.o
lines_patch_mod.o           : ray_patch_mod.o
back_line_patch_mod.o       : ray_patch_mod.o composition_mod.o
back_opac_patch_mod.o       : ray_patch_mod.o composition_mod.o
cont_opac_patch_mod.o       : ray_patch_mod.o atom_mod.o
net_opac_patch_mod.o        : ray_patch_mod.o lines_patch_mod.o back_line_patch_mod.o cont_opac_patch_mod.o
out_rad_patch_mod.o         : ray_patch_mod.o rad_patch_mod.o
proj_patch_mod.o            : ray_patch_mod.o net_opac_patch_mod.o
mean_rad_patch_mod.o        : ray_patch_mod.o rad_patch_mod.o proj_patch_mod.o
mmath.o                     : refsor_mod.o unista_mod.o mconst.o prec_mod.o
