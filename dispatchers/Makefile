VPATH := $(VPATH):$(TOP)/dispatchers

dispatcher_mod.o	: io_mod.o trace_mod.o mpi_mod.o dispatcher0_mod.o dispatcher1_mod.o dispatcher2_mod.o dispatcher3_mod.o
dispatcher0_mod.o	: io_unit_mod.o trace_mod.o task_list_mod.o
dispatcher1_mod.o	: io_unit_mod.o trace_mod.o task_list_mod.o process_mod.o global_mod.o
dispatcher2_mod.o	: io_unit_mod.o trace_mod.o task_list_mod.o
dispatcher3_mod.o	: io_unit_mod.o trace_mod.o task_list_mod.o global_mod.o mpi_mod.o mpi_io_mod.o
