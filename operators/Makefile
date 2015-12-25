# Append to VPATH
VPATH := $(VPATH):$(TOP)/operators

# Dependencies
bits_mod.o              : io_mod.o trace_mod.o
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
tensor_mod.o            : io_mod.o trace_mod.o vector_mod.o
patches_mod.o           : io_mod.o trace_mod.o bits_mod.o
hd_mod.o                : io_mod.o trace_mod.o scalar_mod.o vector_mod.o tensor_mod.o stagger_mod.o
ghd_mod.o               : io_mod.o trace_mod.o hd_mod.o laplace_mod.o stagger_mod.o boundary_mod.o
laplace_mod.o           : io_mod.o trace_mod.o patches_mod.o mpi_mod.o boundary_mod.o
timestep_mod.o          : io_mod.o trace_mod.o $(PDE)_mod.o patches_mod.o
interpolation_mod.o     : io_mod.o trace_mod.o patches_mod.o
patch_list_mod.o        : io_mod.o timestep_mod.o refine_mod.o
refine_mod.o            : io_mod.o trace_mod.o timestep_mod.o
stagger_mod.o           : io_mod.o trace_mod.o
