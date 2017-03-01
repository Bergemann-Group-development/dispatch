# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Common configuration
UPD = adams_bashforth
SOLVER = zeus3d
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

test.x: io_mod.o task_mod.o patch_mod.o volleyball_mod.o link_mod.o
experiment_mod.o: io_mod.o trace_mod.o mhd_mod.o
