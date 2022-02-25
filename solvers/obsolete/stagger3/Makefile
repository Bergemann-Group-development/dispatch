
# Dependencies
mhd_mod.o       : io_mod.o trace_mod.o mpi_mod.o gpatch_mod.o scalar_mod.o vector_mod.o \
                  tensor_mod.o mesh_mod.o initial_mod.o download_mod.o index_mod.o timestep_mod.o
stagger_mod.o	: io_mod.o mesh_mod.o
initial_mod.o   : io_mod.o trace_mod.o mesh_mod.o random_mod.o index_mod.o
scalar_mod.o	: io_mod.o trace_mod.o mesh_mod.o stagger_mod.o mesh_mod.o
vector_mod.o	: io_mod.o trace_mod.o mesh_mod.o stagger_mod.o mesh_mod.o scalar_mod.o
tensor_mod.o	: io_mod.o trace_mod.o mesh_mod.o stagger_mod.o mesh_mod.o
force_mod.o     : io_mod.o trace_mod.o mesh_mod.o vector_mod.o
index_mod.o	: 

