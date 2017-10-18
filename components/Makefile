# Set search path
VPATH := $(VPATH):$(TOP)/components

# Dependencies
cartesian_mod.o         : io_mod.o trace_mod.o task_mod.o task_list_mod.o experiment_mod.o mpi_coords_mod.o boundaries_mod.o
rubiks_cube_mod.o       : io_mod.o trace_mod.o task_mod.o task_list_mod.o experiment_mod.o
rubiks_mod.o            : io_mod.o trace_mod.o task_mod.o task_list_mod.o experiment_mod.o scene_mod.o
scene_mod.o             :
ray_geometry_mod.o      : io_mod.o trace_mod.o
ray_geometry_v2_mod.o   : io_mod.o trace_mod.o
ray_geometry_v3_mod.o   : io_mod.o trace_mod.o
setup_mod.o             : io_mod.o mpi_mod.o scaling_mod.o opacity_mod.o eos_mod.o refine_mod.o download_mod.o
volleyball_mod.o        : io_mod.o trace_mod.o bits_mod.o mesh_mod.o task_mod.o task_list_mod.o experiment_mod.o
