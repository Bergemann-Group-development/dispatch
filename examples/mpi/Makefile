# Make $(PROG).x from $(PROG).f90
default: minmax_example.x coords_example.x send_example.x all_example.x \
         boundary_example.x io_example.x put_example.x get_example.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = \
omp_mod.o \
mpi_mod.o \
mpi_coord_mod.o \
mpi_reduce_mod.o \
mpi_all_mod.o \
mpi_par_io_mod.o \
mpi_rma_mod.o \
mpi_send_mod.o \
io_mod.o \
mpi_timer.o

# Dependencies
minmax_example.x        : $(OBJ)
send_example.x          : $(OBJ)
get_example.x           : $(OBJ)
put_example.x           : $(OBJ)
coords_example.x        : $(OBJ)
all_example.x           : $(OBJ)
boundary_example.x      : $(OBJ)
par_io_example.x            : $(OBJ)
