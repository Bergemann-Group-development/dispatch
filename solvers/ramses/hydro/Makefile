# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)/..
include $(TOP)/solvers/ramses/Makefile

# Dependencies
riemann_mod.o : parameters_mod.o io_mod.o trace_mod.o eos_mod.o
mhd_mod.o     : parameters_mod.o extras_mod.o riemann_mod.o link_mod.o \
                index_mod.o
rmhd_mod.o    : mhd_mod.o rt_mod.o
