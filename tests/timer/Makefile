# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

test.x: omp_mod.o omp_timer_mod.o omp_lock_mod.o io_unit_mod.o
