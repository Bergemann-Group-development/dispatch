# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)

# Dependencies
mhd_rt_mod.o            : mhd_mod.o trace_mod.o microphysics_mod.o rt_data_mod.o rt_internal_mod.o
rt_mod.o                : ray_geometry_v3_mod.o mhd_mod.o link_mod.o rt_boundary_mod.o mesh_mod.o
rt_data_mod.o           : microphysics_mod.o rt_mod.o
rt_internal_mod.o       : microphysics_mod.o rt_mod.o rt_data_mod.o
radiative_transfer_mod.o: ray_geometry_v3_mod.o microphysics_mod.o rt_mod.o list_mod.o \
                          rt_boundary_mod.o rt_data_mod.o rt_internal_mod.o \
                          mhd_mod.o mhd_rt_mod.o experiment_mod.o task_list_mod.o
rt_boundary_mod.o       : io_mod.o trace_mod.o units_mod.o scaling_mod.o
