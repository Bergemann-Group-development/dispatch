# Vpath
VPATH := $(VPATH):$(TOP)/solvers/ramses

# Dependencies
hydro_mod.o      :
index_mod.o      : io_mod.o
initial_mod.o    : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o \
                   index_mod.o units_mod.o parameters_mod.o
force_mod.o      : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o \
                   mpi_global_mod.o
