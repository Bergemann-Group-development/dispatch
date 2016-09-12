# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/rt_$(RT_SOLVER)

# Dependencies
mhd_rt_mod.o            : mhd_mod.o trace_mod.o microphysics_mod.o
rt_mod.o                : ray_geometry_v2_mod.o experiment_mod.o link_mod.o rt_boundary_mod.o
rt_data_mod.o           : microphysics_mod.o rt_mod.o
rt_internal_mod.o       : microphysics_mod.o rt_mod.o rt_data_mod.o
radiative_transfer_mod.o: ray_geometry_v2_mod.o microphysics_mod.o rt_mod.o list_mod.o rt_boundary_mod.o rt_data_mod.o rt_internal_mod.o

