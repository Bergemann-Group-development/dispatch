
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
$(EXEC) : setup_mod.o
setup_mod.o:		list_mod.o task_list_mod.o task_mod.o io_mod.o trace_mod.o experiment_mod.o components_mod.o scaling_mod.o parameters_mod.o rubiks_mod.o
pboundary_mod.o: 	io_mod.o trace_mod.o kinds_mod.o scaling_mod.o mesh_mod.o mhd_mod.o parameters_mod.o boundary_wavekill_mod.o
experiment_mod.o:	io_mod.o trace_mod.o mesh_mod.o solver_mod.o pboundary_mod.o potential_mod.o scaling_mod.o parameters_mod.o index_mod.o
potential_mod.o:	io_mod.o trace_mod.o scaling_mod.o mesh_mod.o parameters_mod.o kinds_mod.o
scaling_mod.o:		io_mod.o trace_mod.o
parameters_mod.o:	io_mod.o trace_mod.o scaling_mod.o
components_mod.o:	io_mod.o trace_mod.o scaling_mod.o task_mod.o task_list_mod.o mhd_mod.o mesh_mod.o bits_mod.o experiment_mod.o parameters_mod.o
initial_mod.o:		io_mod.o trace_mod.o mesh_mod.o kinds_mod.o
force_mod.o:		io_mod.o trace_mod.o scaling_mod.o mesh_mod.o parameters_mod.o kinds_mod.o potential_mod.o stagger_mod.o vector_ops.o
