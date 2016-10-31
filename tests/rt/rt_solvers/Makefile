# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

def: test2

# Common configuration
EXEC = test2.x
TOP = ../../..
include $(TOP)/config/Makefile

<<<<<<< HEAD
test2: test2.x
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test2.x ) | tee -a performance/$(HOST).log
=======
test: test.x
	mkdir -p performance
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test.x ) | tee -a performance/$(HOST).log
>>>>>>> edc2ad4b6ccc3f3c208ed674e587ebda3c5f7528

test2.x: omp_mod.o omp_timer_mod.o ray_geometry_v3_mod.o
