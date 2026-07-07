# Dependencies
read_external_mod.o         : io_mod.o mesh_mod.o task_list_mod.o experiment_mod.o \
                              trace_mod.o task_mod.o link_mod.o bits_mod.o halo_mod.o scaling_mod.o \
                              units_mod.o timer_mod.o patch_mod.o list_mod.o \
                              ray_patch_mod.o read_atmos_mod.o
read_atmos_mod.o            : io_mod.o os_mod.o scaling_mod.o m3d_patch_mod.o
experiment_mod.o            : trace_mod.o io_mod.o kinds_mod.o m3d_halo_mod.o composition_patch_mod.o read_atmos_mod.o
mmath_mod.o                 : refsor_mod.o const_mod.o prec_mod.o search_mod.o
atom_mod.o                  : m3d_debug_mod.o composition_mod.o trans_mod.o op_math_mod.o spectrum_mod.o element_mod.o molecule_mod.o handyfuncs_mod.o hprofl_mod.o
spectrum_mod.o              : io_mod.o timer_mod.o trace_mod.o const_mod.o trans_mod.o handyfuncs_mod.o
m3d_halo_mod.o              : halo_mod.o gpatch_mod.o handyfuncs_mod.o
m3d_patch_mod.o             : gpatch_mod.o atom_mod.o ray_mod.o m3d_halo_mod.o handyfuncs_mod.o m3d_debug_mod.o composition_mod.o
ray_patch_mod.o             : m3d_patch_mod.o
rt_patch_mod.o              : composition_mod.o ray_patch_mod.o
nlte_patch_mod.o            : rt_patch_mod.o
composition_patch_mod.o     : nlte_patch_mod.o rt_patch_mod.o hvcs_mod.o
bftrans_mod.o               : os_mod.o
bbtrans_mod.o               : handyfuncs_mod.o trace_mod.o
ray_mod.o                   : mmath_mod.o radau_mod.o
prec_mod.o                  : const_mod.o os_mod.o
patch_mod.o                 : os_mod.o handyfuncs_mod.o
mpi_debug_mod.o             : mpi_mod.o io_mod.o
setup_mod.o                 : io_mod.o trace_mod.o mpi_mod.o omp_mod.o omp_lock_mod.o \
                              mpi_mesg_mod.o timer_mod.o scaling_mod.o
boundaries_mod.o            : io_mod.o trace_mod.o mpi_mod.o bits_mod.o mesh_mod.o \
                              index_mod.o patch_mod.o
gpatch_mod.o                : boundaries_mod.o
m3d_debug_mod.o             :
refsor_mod.o                :
search_mod.o                :
dispatch.o                  : io_mod.o 