# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

def: test2

# Common configuration
EXEC = test2.x
TOP = ../../..
include $(TOP)/config/Makefile

test2: test2.x
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test2.x ) | tee -a performance/$(HOST).log

test2.x: omp_mod.o omp_timer_mod.o ray_geometry_v3_mod.o
