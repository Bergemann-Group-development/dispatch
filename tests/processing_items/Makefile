# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Optional extra object files
dispatch.x: $(OBJ)
test.x: $(OBJ)
OBJ = bits_mod.o item_mod.o task_mod.o list_mod.o task_list_mod.o experiment_mod.o

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: task_list_mod.o

dispatch.x: experiment_mod.o
experiment_mod.o: task_list_mod.o
