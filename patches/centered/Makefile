# Append to VPATH
VPATH := $(VPATH):$(TOP)/patches/$(MESH)

# Local object files to include (hardwired later)
OBJ += \
$(PROBLEM)_mod.o

# Dependencies
$(PROBLEM)_mod.o 	: patch_mod.o mpi_mod.o
