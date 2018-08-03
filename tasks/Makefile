# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Dependencies
task_mod.o:             io_mod.o mesh_mod.o bits_mod.o mpi_mod.o omp_lock_mod.o scene_mod.o mpi_mesg_mod.o random_mod.o dll_mod.o
patch_mod.o:            task_mod.o mesh_mod.o mpi_mesg_mod.o link_mod.o omp_timer_mod.o \
                        timer_mod.o boundaries_mod.o
gpatch_mod.o:           io_mod.o trace_mod.o patch_mod.o download_mod.o initial_mod.o kinds_mod.o data_hub_mod.o \
                        poisson_mod.o mesh_mod.o scaling_mod.o data_io_mod.o units_mod.o
refine_mod.o:           io_mod.o trace_mod.o bits_mod.o list_mod.o link_mod.o kinds_mod.o scaling_mod.o \
                        download_mod.o mesh_mod.o solver_mod.o patch_mod.o sinkparticles_mod.o
random_mod.o:           io_mod.o trace_mod.o
radiation_mod.o:        io_mod.o trace_mod.o mesh_mod.o link_mod.o
sinkparticles_mod.o:	io_unit_mod.o trace_mod.o download_mod.o data_hub_mod.o link_mod.o list_mod.o particle_history_mod.o
index_mod.o:            io_mod.o
