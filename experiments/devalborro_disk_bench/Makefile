
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Path to ZEUScore shared library.
# Don't forget to update the LD_LIBRARY_PATH environment variable.
LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore

# Targets and methods
UPD = adams_bashforth
SOLVER = zeus3d
EXEC = dispatch.x
TOP = ../..
sinclude $(TOP)/config/Makefile

# Dependencies
$(EXEC) : setup_mod.o dispatcher_mod.o download_mod.o scaling_mod.o mpi_mod.o io_mod.o
setup_mod.o:		list_mod.o task_list_mod.o task_mod.o io_mod.o trace_mod.o experiment_mod.o components_mod.o scaling_mod.o parameters_mod.o
pboundary_mod.o: 	io_mod.o trace_mod.o kinds_mod.o scaling_mod.o mesh_mod.o solver_mod.o parameters_mod.o boundary_wavekill_mod.o
experiment_mod.o:	io_mod.o trace_mod.o mesh_mod.o solver_mod.o pboundary_mod.o potential_mod.o scaling_mod.o parameters_mod.o
potential_mod.o:	io_mod.o trace_mod.o scaling_mod.o mesh_mod.o parameters_mod.o solver_mod.o kinds_mod.o
scaling_mod.o:		io_mod.o trace_mod.o
parameters_mod.o:	io_mod.o trace_mod.o scaling_mod.o
components_mod.o:	io_mod.o trace_mod.o scaling_mod.o task_mod.o task_list_mod.o solver_mod.o mesh_mod.o bits_mod.o experiment_mod.o parameters_mod.o
initial_mod.o:		io_mod.o trace_mod.o mesh_mod.o kinds_mod.o index_mod.o
