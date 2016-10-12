# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

test: test.x
	mkdir -p performance
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test.x ) | tee -a performance/$(HOST).log

poisson_mod.o   : io_unit_mod.o trace_mod.o mesh_mod.o
test.x          : io_mod.o io_unit_mod.o trace_mod.o mesh_mod.o patch_mod.o poisson_mod.o
