# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

cg_mod.o        : lagrange_mod.o timer_mod.o
test.x          : cg_mod.o
