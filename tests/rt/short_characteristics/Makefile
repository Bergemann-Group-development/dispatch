# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

$(EXEC)		: io_mod.o mpi_mod.o rt_sc_mod.o
rt_sc_mod.o	: io_mod.o trace_mod.o timer_mod.o
