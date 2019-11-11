# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Dependencies
task_mod.o:       io_mod.o mesh_mod.o bits_mod.o mpi_mod.o omp_lock_mod.o \
                  mpi_mesg_mod.o random_mod.o dll_mod.o aux_mod.o timer_mod.o
patch_mod.o:      task_mod.o mpi_mesg_mod.o link_mod.o omp_timer_mod.o \
                  timer_mod.o boundaries_mod.o connect_mod.o omp_mod.o \
                  shared_mod.o
gpatch_mod.o:     patch_mod.o initial_mod.o kinds_mod.o mesh_mod.o scaling_mod.o \
                  data_io_mod.o units_mod.o guard_zones_mod.o download_mod.o \
                  data_hub_mod.o validate_mod.o list_mod.o
refine_mod.o:     io_mod.o trace_mod.o bits_mod.o list_mod.o link_mod.o kinds_mod.o scaling_mod.o \
                  download_mod.o mesh_mod.o solver_mod.o patch_mod.o math_mod.o \
                  omp_lock_mod.o data_io_mod.o shared_mod.o
random_mod.o:     io_mod.o trace_mod.o
radiation_mod.o:  io_mod.o trace_mod.o mesh_mod.o link_mod.o
index_mod.o:      io_mod.o trace_mod.o
cpatch_mod.o:     io_mod.o trace_mod.o solver_mod.o
