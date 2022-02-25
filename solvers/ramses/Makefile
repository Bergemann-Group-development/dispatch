# Vpath
VPATH := $(VPATH):$(TOP)/solvers/ramses

# Dependencies
parameters_mod.o :
slope_mod.o      : patch_mod.o parameters_mod.o
solver_mod.o     : mhd_mod.o kinds_mod.o validate_mod.o stagger_2nd.o \
                   parameters_mod.o
initial_mod.o    : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o \
                   index_mod.o units_mod.o parameters_mod.o eos_mod.o
force_mod.o      : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o \
                   mpi_global_mod.o
