VPATH := $(VPATH):$(TOP)/dispatchers

dispatcher_mod.o        : io_mod.o trace_mod.o mpi_mod.o list_mod.o task_list_mod.o \
                          dispatcher0_mod.o dispatcher1_mod.o dispatcher2_mod.o \
                          dispatcher3_mod.o dispatcher4_mod.o dispatcher5_mod.o \
                          dispatcher6_mod.o experiment_mod.o
dispatcher0_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o data_io_mod.o omp_lock_mod.o validate_mod.o
dispatcher1_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o process_mod.o global_mod.o
dispatcher2_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o
dispatcher3_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o global_mod.o mpi_mod.o mpi_io_mod.o
dispatcher4_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o global_mod.o mpi_mod.o mpi_io_mod.o
dispatcher5_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o omp_lock_mod.o
dispatcher6_mod.o       : io_unit_mod.o trace_mod.o task_list_mod.o omp_lock_mod.o
