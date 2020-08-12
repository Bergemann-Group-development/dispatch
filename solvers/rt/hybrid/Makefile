VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)

rt_mod.o                 : io_mod.o trace_mod.o index_mod.o link_mod.o \
                           eos_mod.o gpatch_mod.o math_mod.o bits_mod.o rt_solver_mod.o \
                           scalar_mod.o timestep_mod.o rt_boundaries_mod.o rt_components_mod.o \
                           scratch_mod.o
rt_boundaries_mod.o      : io_mod.o trace_mod.o kinds_mod.o mpi_mod.o bits_mod.o \
                           mesh_mod.o eos_mod.o units_mod.o scaling_mod.o ray_geometry_v3_mod.o
rt_solver_mod.o          : io_mod.o trace_mod.o kinds_mod.o
rt_components_mod.o      : patch_mod.o ray_geometry_v3_mod.o trace_mod.o io_mod.o link_mod.o

