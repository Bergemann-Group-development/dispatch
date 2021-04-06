# The OFFLOAF macro determines which device_mod.o to use

OFFLOAD ?= CPU

VPATH := $(VPATH):$(TOP)/offload:$(TOP)/offload/$(OFFLOAD)
sinclude $(TOP)/offload/$(OFFLOAD)/Makefile

offload_mod.o :  io_mod.o trace_mod.o list_mod.o bits_mod.o patch_mod.o \
                 omp_timer_mod.o
