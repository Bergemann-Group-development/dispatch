# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

def: test

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

test: test.x
	mkdir -p performance
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test.x ) | tee -a performance/$(HOST).log

test.x: omp_mod.o omp_timer_mod.o random_mod.o
