# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

OBJ += timer_mod.o
default: test_defaults.x test_reduction.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test_defaults.x: $(OBJ)
test_reduction.x: $(OBJ)
