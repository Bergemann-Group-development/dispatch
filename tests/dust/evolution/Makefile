# Option groups
OPTS = full_debug
#OPTS = debug
#OPTS =
#OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

LIBS+=-llapack 

test.x: units_mod.o math_mod.o
