# Make $(PROG).x from $(PROG).f90
default: minmax_example.x coords_example.x send_example.x all_example.x \
         boundary_example.x io_example.x put_example.x get_example.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = omp_mod.o mpi_basic.o mpi_coords.o mpi_reduce.o mpi_all.o mpi_par_io.o mpi_rma.o mpi_send.o io_mod.o mpi_timer.o

# Dependencies
minmax_example.x        : $(OBJ)
send_example.x          : $(OBJ)
get_example.x           : $(OBJ)
put_example.x           : $(OBJ)
coords_example.x        : $(OBJ)
all_example.x           : $(OBJ)
boundary_example.x      : $(OBJ)
io_example.x            : $(OBJ)
