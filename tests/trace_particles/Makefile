# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

OBJ += io_mod.o omp_timer_mod.o dll_mod.o random_mod.o
default: test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: $(OBJ)
