# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
EXEC = test.x
UPD = adams_bashforth
SOLVER = stagger2
TOP = ../../..
include $(TOP)/config/Makefile

$(EXEC)		: io_mod.o mpi_mod.o rt_sc_mod.o rt_patch_mod.o
rt_sc_mod.o	: io_mod.o trace_mod.o timer_mod.o
rt_patch_mod.o  : io_mod.o trace_mod.o radau_mod.o mhd_mod.o index_mod.o
