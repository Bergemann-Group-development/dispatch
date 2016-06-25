# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

test_queue_by_time.o: omp_lock_mod.o

test.x: omp_lock_mod.o random_mod.o
