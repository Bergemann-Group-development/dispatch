# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

OBJ += list_mod.o task_mod.o patch_mod.o mhd_mod.o timer_mod.o
default: test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: $(OBJ)
