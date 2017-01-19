
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Targets and methods
UPD = adams_bashforth
SOLVER = zeus3d
EXEC = dispatch.x
TOP = ../..
sinclude $(TOP)/config/Makefile

# Dependencies
$(EXEC) : setup_mod.o
setup_mod.o:		task_list_mod.o io_mod.o trace_mod.o experiment_mod.o components_mod.o scaling_mod.o parameters_mod.o rubiks_mod.o timestep_mod.o
scaling_mod.o:		io_mod.o trace_mod.o cgs_mod.o
planetesimal_mod.o:	io_mod.o trace_mod.o patch_mod.o mhd_mod.o
experiment_mod.o:	io_mod.o trace_mod.o timestep_mod.o mesh_mod.o mhd_mod.o boundary_mod.o force_mod.o
components_mod.o:	io_mod.o trace_mod.o task_mod.o mhd_mod.o task_list_mod.o experiment_mod.o bits_mod.o parameters_mod.o
boundary_mod.o:		io_mod.o trace_mod.o scaling_mod.o mesh_mod.o mhd_mod.o parameters_mod.o
force_mod.o:		io_mod.o trace_mod.o scaling_mod.o mesh_mod.o mhd_mod.o parameters_mod.o kinds_mod.o
parameters_mod.o:	io_mod.o trace_mod.o scaling_mod.o cgs_mod.o planetesimal_mod.o
