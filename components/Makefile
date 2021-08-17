# Set search path
VPATH := $(VPATH):$(TOP)/components
sinclude $(TOP)/components/external/Makefile

# Dependencies
cartesian_mod.o         : io_mod.o trace_mod.o experiment_mod.o  task_list_mod.o \
                          mpi_coords_mod.o omp_mod.o
mesh_reader.o           : io_mod.o trace_mod.o
rubiks_cube_mod.o       : io_mod.o trace_mod.o task_mod.o task_list_mod.o experiment_mod.o
rubiks_mod.o            : io_mod.o trace_mod.o task_mod.o task_list_mod.o experiment_mod.o shared_mod.o
shared_mod.o            :
ray_geometry_mod.o      : io_mod.o trace_mod.o
ray_geometry_v2_mod.o   : io_mod.o trace_mod.o
ray_geometry_v3_mod.o   : io_mod.o trace_mod.o
read_snapshot_mod.o     : task_list_mod.o lagrange_mod.o experiment_mod.o bifrost_scaling_mod.o scaling_mod.o
setup_mod.o             : io_mod.o trace_mod.o mpi_mod.o omp_mod.o omp_lock_mod.o \
                          mpi_mesg_mod.o timer_mod.o eos_mod.o opacity_mod.o scaling_mod.o \
			  download_mod.o refine_mod.o
volleyball_mod.o        : io_mod.o trace_mod.o bits_mod.o mesh_mod.o task_mod.o task_list_mod.o experiment_mod.o
process_mod.o           : trace_mod.o link_mod.o task_list_mod.o mpi_comm_mod.o
ramses_mod.o            : io_mod.o trace_mod.o scaling_mod.o mhd_mod.o force_mod.o initial_mod.o units_mod.o \
                          experiment_mod.o star_mod.o task_list_mod.o index_mod.o dll_mod.o \
			  hash_table_mod.o omp_timer_mod.o shared_mod.o
cylinder_mod.o          : io_mod.o io_unit_mod.o trace_mod.o task_list_mod.o solver_mod.o
sphere_points_mod.o     : io_mod.o io_unit_mod.o trace_mod.o
volleyball_layers_mod.o : sphere_points_mod.o dll_mod.o
cartesian_vb_mod.o      : io_mod.o trace_mod.o experiment_mod.o  task_list_mod.o \
                          mpi_coords_mod.o omp_mod.o volleyball_layers_mod.o
restart_mod.o           : io_mod.o trace_mod.o task_list_mod.o task_factory_mod.o \
                          amr_io_mod.o read_stagger_mod.o read_bifrost_mod.o
read_stagger_mod.o      : task_list_mod.o lagrange_mod.o experiment_mod.o scaling_mod.o
read_bifrost_mod.o      : task_list_mod.o lagrange_mod.o experiment_mod.o scaling_mod.o \
                          bifrost_scaling_mod.o
