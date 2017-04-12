# Set search path
VPATH := $(VPATH):$(TOP)/tasks

# Dependencies
item_mod.o:             bits_mod.o trace_mod.o
task_mod.o:             io_mod.o mesh_mod.o bits_mod.o mpi_mod.o omp_lock_mod.o scene_mod.o mpi_mesg_mod.o
patch_mod.o:            task_mod.o mesh_mod.o buffered_output_mod.o mpi_mesg_mod.o link_mod.o omp_timer_mod.o timer_mod.o boundaries_mod.o
gpatch_mod.o:           io_mod.o trace_mod.o patch_mod.o download_mod.o initial_mod.o kinds_mod.o
hd_mod.o:               patch_mod.o trace_mod.o
unigrid_mod.o:          patch_mod.o mpi_coords_mod.o
refine_mod.o:           io_mod.o trace_mod.o patch_mod.o experiment_mod.o
scaling_mod.o:          io_mod.o cgs_mod.o
download_mod.o:         io_mod.o trace_mod.o task_mod.o patch_mod.o link_mod.o interpolate_mod.o kinds_mod.o
random_mod.o:           io_mod.o trace_mod.o
lagrange_mod.o:         io_mod.o
radiation_mod.o:        io_mod.o trace_mod.o mesh_mod.o link_mod.o
cgs_mod.o:
interpolate_mod.o:      io_mod.o trace_mod.o mpi_mod.o kinds_mod.o
