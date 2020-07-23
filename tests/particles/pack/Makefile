# Option groups
#OPTS = full_debug
#OPTS = debug
#OPTS =
OPTS = optimized

EXEC = test.x
default: $(EXEC)

# Common configuration
TOP = ../../..
include $(TOP)/config/Makefile

$(EXEC): io_mod.o trace_mod.o mpi_buffer_mod.o particles_mod.o
particles_mod.o: io_mod.o trace_mod.o random_mod.o
