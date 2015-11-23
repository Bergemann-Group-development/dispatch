# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Local object files to include (hardwired later)
ifeq ($(EXEC),dispatch.x)
  OBJ += \
  io_mod.o \
  trace_mod.o
endif

# Dependencies
trace_mod.o:            io_mod.o mpi_timer_mod.o
io_mod.o:               mpi_mod.o mpi_timer_mod.o

