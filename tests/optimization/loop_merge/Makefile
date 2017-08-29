# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
SOLVER = stagger2
STAGGER = stagger_4th.o
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

test.x: $(STAGGER) omp_timer_mod.o
