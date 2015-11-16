# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target
default: run

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = mpi_timer.o
classical.x: $(OBJ)
modules.x: $(OBJ)
objects.x: $(OBJ)

# Dependencies
run: classical.x modules.x objects.x
	./classical.x  < input.nml | tee classical.log
	./modules.x    < input.nml | tee modules.log
	./objects.x    < input.nml | tee objects.log
