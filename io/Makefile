# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Local object files to include (hardwired later)
#OBJ += \
#mesh_mod.o

# Dependencies
trace_mod.o:            io_mod.o mpi_timer.o
io_mod.o:               mpi_basic.o mpi_timer.o
