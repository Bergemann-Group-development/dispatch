
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
$(EXEC) : setup_mod.o dispatcher_mod.o experiment_setup_mod.o
experiment_mod.o:	io_mod.o trace_mod.o solver_mod.o pboundaries_mod.o pebble_patch_mod.o \
					potential_mod.o kinds_mod.o
scaling_mod.o:		io_mod.o trace_mod.o
experiment_setup_mod.o:	io_mod.o trace_mod.o scaling_mod.o task_list_mod.o mesh_mod.o bits_mod.o \
						experiment_mod.o math_mod.o mpi_mod.o task_mod.o bits_mod.o mesh_mod.o
initial_mod.o:		io_mod.o trace_mod.o mesh_mod.o kinds_mod.o index_mod.o math_mod.o units_mod.o
pboundaries_mod.o:	io_mod.o trace_mod.o mesh_mod.o bits_mod.o index_mod.o kinds_mod.o \
					mpi_mod.o patch_mod.o vector_ops.o
potential_mod.o:	io_mod.o trace_mod.o mpi_mod.o mesh_mod.o kinds_mod.o index_mod.o scaling_mod.o
pebble_mod.o:
pebbles_mod.o:		io_mod.o trace_mod.o bits_mod.o pebble_mod.o omp_lock_mod.o
pebble_trace_mod.o:	io_mod.o trace_mod.o bits_mod.o omp_mod.o omp_lock_mod.o pebble_mod.o
pebble_patch_mod.o:	io_mod.o trace_mod.o bits_mod.o omp_mod.o omp_lock_mod.o omp_timer_mod.o link_mod.o \
					kinds_mod.o task_mod.o patch_mod.o scaling_mod.o units_mod.o pebble_mod.o \
					pebbles_mod.o pebble_trace_mod.o mesh_mod.o solver_mod.o index_mod.o random_mod.o
