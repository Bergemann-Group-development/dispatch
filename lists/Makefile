# Append to VPATH
VPATH := $(VPATH):$(TOP)/lists

# Dependencies (should only be listed for files that exist in this directory)
sll_mod.o           :
dll_mod.o           : io_mod.o trace_mod.o omp_lock_mod.o
counters_mod.o      : io_mod.o trace_mod.o omp_lock_mod.o
plist_mod.o         : dll_mod.o omp_timer_mod.o
link_mod.o          : task_mod.o mpi_mesg_mod.o io_unit_mod.o timer_mod.o \
                      omp_lock_mod.o
task_block_mod.o    : link_mod.o
boundary_mod.o      : link_mod.o
hash_table_mod.o    : link_mod.o shared_mod.o os_mod.o
guards_mod.o        : link_mod.o patch_mod.o lagrange_mod.o interp_mod.o conserve_mod.o \
                      math_mod.o boundary_mod.o validate_mod.o
halo_mod.o          : guards_mod.o boundary_mod.o region_mod.o
list_mod.o          : halo_mod.o patch_mod.o hash_table_mod.o running_mean_mod.o \
                      mpi_window_mod.o os_mod.o
mulligan_mod.o      : list_mod.o
load_balance_mod.o  : list_mod.o experiment_mod.o
task_mesg_mod.o     : load_balance_mod.o os_mod.o mpi_lib.o
task_list_mod.o     : task_mesg_mod.o refine_mod.o mulligan_mod.o offload_mod.o \
                      memory_mod.o halo_mod.o os_mod.o

