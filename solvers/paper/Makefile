# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/paper

PDE     ?= mhd_mod.o
STAGGER ?= stagger_mod.o

# Dependencies
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
vector_ops.o            : $(STAGGER) mesh_mod.o
stagger_mod.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o \
                          stagger_16.o stagger_20.o stagger_24.o stagger_32.o stagger_36.o
stagger_4th.o           : io_mod.o trace_mod.o scalar_mod.o vector_mod.o
mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o gpatch_mod.o scaling_mod.o units_mod.o download_mod.o \
                          $(STAGGER) eos_mod.o index_mod.o timestep_mod.o vector_ops.o
hd_mod.o                : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o gpatch_mod.o scaling_mod.o units_mod.o download_mod.o \
                          $(STAGGER) eos_mod.o index_mod.o timestep_mod.o vector_ops.o
patches_mod.o           : io_mod.o trace_mod.o timestep_mod.o
initial_mod.o           : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o $(STAGGER) index_mod.o
force_mod.o             : io_mod.o trace_mod.o random_mod.o mesh_mod.o vector_mod.o $(STAGGER) scaling_mod.o
boundary_mod.o          : io_mod.o trace_mod.o mesh_mod.o
mhd_rt_mod.o            : $(PDE) trace_mod.o microphysics_mod.o
gmhd_mod.o              : io_mod.o omp_timer_mod.o trace_mod.o task_mod.o patch_mod.o scalar_mod.o laplace_mod.o $(PDE) boundary_mod.o scalar_mod.o
rmhd_mod.o              : io_mod.o trace_mod.o $(PDE) rt_mod.o
