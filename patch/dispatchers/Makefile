VPATH := $(VPATH):$(TOP)/dispatchers

dispatcher_mod.o		: io_mod.o mpi_mesg_mod.o dispatcher_method1_mod.o dispatcher_method2_mod.o
dispatcher_method1_mod.o	: io_unit_mod.o trace_mod.o link_mod.o task_list_mod.o
dispatcher_method2_mod.o	: io_unit_mod.o trace_mod.o link_mod.o task_list_mod.o
none_mod.o			: task_list_mod.o
