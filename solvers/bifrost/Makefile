
# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/bifrost:$(TOP)/solvers/bifrost/shadow

# Connect to the make system of the code
BIFROST = ../../../Bifrost

include $(TOP)/solvers/bifrost/modules.mkf
sinclude modules.mkf
include $(TOP)/solvers/bifrost/bifrost.mkf

# Dependencies
initial_mod.o		: io_mod.o trace_mod.o mesh_mod.o index_mod.o
extern_mod.o		: params.o mpi.o timing_empty.o main.o solver_mod.o
solver_mod.o		: io_mod.o patch_mod.o download_mod.o bifrost_mod.o mesh_interface_mod.o bifrost_rt_mod.o
bifrost_mod.o       : main.o params.o mpi.o $(UNITS).o stagger_mesh_$(MP).o $(TIMESTEP).o index_mod.o mesh_interface_mod.o \
                      mesh_mod.o
mesh_interface_mod.o	: trace_mod.o params.o stagger_mesh_$(MP).o io_mod.o io_unit_mod.o mesh_mod.o
bifrost_rt_mod.o	: trace_mod.o params.o $(RADIATION).o
mpi.o			: trace_mod.o
genrad_$(MP).o		: incrad_$(MP).o
bc_lower_magnetic.o	: spitzer_conductivity_$(MP).o
