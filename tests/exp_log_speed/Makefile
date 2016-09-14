# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized_diagnostics

# Common configuration
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

cases.o: omp_mod.o omp_timer_mod.o
test.x: cases.o

perflog: test.x
	mkdir -p performance
	uname -a | tee -a performance/$(HOST).log
	./test.x | tee -a performance/$(HOST).log

