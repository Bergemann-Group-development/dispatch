# Option groups
OPTS = full_debug
#OPTS = debug
#OPTS = optimized

# Common configuration
SOLVER = ramses/hydro
EOS = itable_mod.o
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

test.x: omp_mod.o code.o
