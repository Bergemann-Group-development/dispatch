# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
bits_mod.o              : io_mod.o trace_mod.o
scalar_mod.o            : io_mod.o trace_mod.o
vector_mod.o            : io_mod.o trace_mod.o scalar_mod.o
tensor_mod.o            : io_mod.o trace_mod.o vector_mod.o
stagger_mod.o           : io_mod.o trace_mod.o
patches_mod.o           : io_mod.o trace_mod.o bits_mod.o timestep_mod.o
hd_solver_mod.o         : io_mod.o trace_mod.o scalar_mod.o vector_mod.o tensor_mod.o stagger_mod.o
ghd_solver_mod.o        : io_mod.o trace_mod.o hd_solver_mod.o laplace_mod.o stagger_mod.o boundary_mod.o
laplace_mod.o           : io_mod.o trace_mod.o patches_mod.o mpi_mod.o boundary_mod.o
interpolation_mod.o     : io_mod.o trace_mod.o patches_mod.o
refine_mod.o            : io_mod.o trace_mod.o timestep_mod.o
random_mod.o		: io_mod.o
scaling_mod.o           : io_mod.o units_mod.o
