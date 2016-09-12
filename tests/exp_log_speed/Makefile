# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized_diagnostics

# Common configuration
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

test.x: omp_mod.o omp_timer_mod.o

perflog: test.x
	mkdir -p performance
	uname -a | tee -a performance/$(HOST).log
	./test.x | tee -a performance/$(HOST).log
