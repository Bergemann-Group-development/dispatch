# Append to VPATH
VPATH := $(VPATH):$(TOP)/lists

# Dependencies
dispatcher_mod.o    : omp_mod.o timer_mod.o io_mod.o io_unit_mod.o trace_mod.o link_mod.o list_mod.o task_list_mod.o mpi_mesg_mod.o
load_balance_mod.o  : io_mod.o trace_mod.o mpi_mesg_mod.o omp_timer_mod.o patch_mod.o random_mod.o omp_lock_mod.o \
                      list_mod.o
link_mod.o          : io_mod.o trace_mod.o task_mod.o mpi_mesg_mod.o io_unit_mod.o timer_mod.o
list_mod.o          : io_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o omp_lock_mod.o \
                      download_mod.o counters_mod.o
task_mesg_mod.o     : omp_timer_mod.o mpi_mod.o list_mod.o task_mod.o bits_mod.o trace_mod.o io_mod.o load_balance_mod.o \
                      experiment_mod.o refine_mod.o omp_lock_mod.o mpi_mesg_mod.o timer_mod.o \
                      direct_io_mod.o io_unit_mod.o hash_table_mod.o
task_list_mod.o     : omp_timer_mod.o mpi_mod.o list_mod.o task_mod.o bits_mod.o \
                      trace_mod.o io_mod.o load_balance_mod.o experiment_mod.o \
                      refine_mod.o omp_lock_mod.o mpi_mesg_mod.o timer_mod.o \
                      direct_io_mod.o data_io_mod.o task_mesg_mod.o io_unit_mod.o \
                      validate_mod.o download_mod.o shared_mod.o
mesg_list_mod.o     : task_list_mod.o
counters_mod.o      : io_mod.o trace_mod.o omp_lock_mod.o
plist_mod.o         : omp_timer_mod.o dll_mod.o
hash_table_mod.o    : mpi_mod.o io_mod.o trace_mod.o
