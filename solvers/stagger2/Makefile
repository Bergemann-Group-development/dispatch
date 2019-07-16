# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/stagger2

PDE     ?= mhd_mod.o
STAGGER ?= stagger_2nd.o

# Dependencies
mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o extras_mod.o scaling_mod.o units_mod.o download_mod.o \
                          $(STAGGER) eos_mod.o index_mod.o timestep_mod.o vector_ops.o kinds_mod.o non_ideal_mod.o \
                          aux_mod.o
hd_mod.o                : io_mod.o trace_mod.o scalar_mod.o vector_mod.o extras_mod.o scaling_mod.o units_mod.o download_mod.o \
                          $(STAGGER) eos_mod.o index_mod.o timestep_mod.o vector_ops.o
patches_mod.o           : io_mod.o trace_mod.o timestep_mod.o
initial_mod.o           : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o $(STAGGER) index_mod.o
force_mod.o             : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o $(STAGGER) scaling_mod.o
boundary_mod.o          : io_mod.o trace_mod.o mesh_mod.o
mhd_rt_mod.o            : $(PDE) trace_mod.o microphysics_mod.o
gmhd_mod.o              : io_mod.o omp_timer_mod.o trace_mod.o task_mod.o patch_mod.o scalar_mod.o laplace_mod.o $(PDE) boundary_mod.o scalar_mod.o
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o eos_mod.o rt_mod.o solver_mod.o timestep_mod.o
solver_mod.o            : io_mod.o trace_mod.o mhd_mod.o eos_mod.o validate_mod.o
non_ideal_mod.o         : io_mod.o trace_mod.o scalar_mod.o $(STAGGER)
index_mod.o		: io_unit_mod.o
