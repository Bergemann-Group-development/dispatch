
# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Targets and methods
UPD = adams_bashforth
SOLVER = stagger2/multifluid
EXEC = dispatch.x
TOP = ../..
sinclude $(TOP)/config/Makefile

# Dependencies
$(EXEC) : setup_mod.o dispatcher_mod.o
setup_mod.o:		list_mod.o task_list_mod.o task_mod.o io_mod.o trace_mod.o experiment_mod.o components_mod.o scaling_mod.o parameters_mod.o rubiks_mod.o
pboundary_mod.o: 	io_mod.o trace_mod.o kinds_mod.o scaling_mod.o units_mod.o mesh_mod.o patch_mod.o parameters_mod.o boundary_wavekill_mod.o
experiment_mod.o:	io_mod.o trace_mod.o mesh_mod.o solver_mod.o pboundary_mod.o potential_mod.o scaling_mod.o parameters_mod.o index_mod.o
potential_mod.o:	io_mod.o trace_mod.o scaling_mod.o mesh_mod.o parameters_mod.o kinds_mod.o
scaling_mod.o:		io_mod.o trace_mod.o parameters_mod.o
parameters_mod.o:	io_mod.o trace_mod.o
components_mod.o:	io_mod.o trace_mod.o scaling_mod.o task_mod.o task_list_mod.o mhd_mod.o mesh_mod.o bits_mod.o experiment_mod.o parameters_mod.o
initial_mod.o:		io_mod.o trace_mod.o mesh_mod.o kinds_mod.o index_mod.o scaling_mod.o units_mod.o parameters_mod.o
force_mod.o:		io_mod.o trace_mod.o scaling_mod.o mesh_mod.o parameters_mod.o kinds_mod.o potential_mod.o stagger_mod.o vector_ops.o
gpatch_mod.o:           io_mod.o trace_mod.o patch_mod.o download_mod.o initial_mod.o kinds_mod.o data_hub_mod.o scalar_mod.o \
                        poisson_mod.o mesh_mod.o units_mod.o scaling_mod.o data_io_mod.o timestep_mod.o pboundary_mod.o os_mod.o

