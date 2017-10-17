# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/rt/hybrid

# Dependencies
rt_solver_mod.o           : io_mod.o trace_mod.o
rt_components_mod.o       : ray_geometry_v3_mod.o trace_mod.o io_mod.o link_mod.o rt_boundaries_mod.o patch_mod.o
rt_mod.o                  : ray_geometry_v3_mod.o patch_mod.o rt_components_mod.o trace_mod.o rt_solver_mod.o eos_mod.o \
                            download_mod.o data_hub_mod.o
radiative_transfer_mod.o  : experiment_mod.o mesh_mod.o ray_geometry_v3_mod.o rt_mod.o
