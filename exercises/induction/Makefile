# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target
default: induction.x

# Common configuration
TOP = ../..
sinclude $(TOP)/config/Makefile
TOP = ../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../../..
sinclude $(TOP)/config/Makefile

OBJ = io_mod.o \
trace_mod.o \
omp_mod.o \
mpi_mod.o \
mpi_timer_mod.o \
mpi_coords_mod.o \
mpi_send_mod.o \
mpi_reduce_mod.o \
mpi_par_io_mod.o \
task_mod.o \
upwind_mod.o \
central_mod.o \
derivatives_mod.o

induction.x: $(OBJ)

# Dependencies
upwind_mod.o: mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o
central_mod.o: mpi_coords_mod.o mpi_send_mod.o derivatives_mod.o
run: induction.x
	./induction.x | tee induction.log

