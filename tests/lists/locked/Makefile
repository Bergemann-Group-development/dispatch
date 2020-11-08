# Option groups
#OPTS = full_debug
OPTS = debug
#OPTS = 
#OPTS = optimized

EXEC = test.x

# Common configuration
TOP = ../../..
include $(TOP)/config/Makefile

dll_mod.o        : random_mod.o trace_mod.o io_mod.o omp_lock_mod.o omp_timer_mod.o
link_mod.o       : dll_mod.o trace_mod.o task_mod.o
list_mod.o       : link_mod.o
task_mesg_mod.o  : list_mod.o
task_list_mod.o  : task_mesg_mod.o

$(EXEC): list_mod.o random_mod.o
