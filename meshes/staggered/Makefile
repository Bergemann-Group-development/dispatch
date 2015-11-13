# Append to VPATH
VPATH := $(VPATH):$(TOP)/meshes/$(MESH)

# Local object files to include (hardwired later)
ifeq ($(EXEC),dispatch.x)
OBJ += \
mesh_mod.o
endif

# Dependencies
mesh_mod.o:             trace_mod.o
