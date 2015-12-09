# Append to VPATH
VPATH := $(VPATH):$(TOP)/operators

# Local object files to include (hardwired later)
ifeq ($(EXEC),dispatch.x)
  OBJ += \
  upwind_mod.o \
  central_mod.o \
  scalar_mod.o \
  vector_mod.o \
  tensor_mod.o \
  derivatives_mod.o
endif

# Dependencies
#central_mod.o           : derivatives_mod.o mpi_coords_mod.o mpi_send_mod.o io_mod.o
#upwind_mod.o            : derivatives_mod.o central_mod.o
#derivatives_mod.o       : io_mod.o
#scalar_mod.o            : derivatives_mod.o trace_mod.o io_mod.o
#vector_mod.o            : scalar_mod.o trace_mod.o io_mod.o
#tensor_mod.o            : vector_mod.o trace_mod.o
scalar_mod.o            : trace_mod.o io_mod.o
vector_mod.o            : scalar_mod.o trace_mod.o io_mod.o
tensor_mod.o            : vector_mod.o trace_mod.o
