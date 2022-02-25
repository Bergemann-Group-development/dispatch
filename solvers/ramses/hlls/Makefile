# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)
include $(TOP)/solvers/ramses/Makefile

# Dependencies
riemann_mod.o : eos_mod.o parameters_mod.o io_mod.o trace_mod.o timer_mod.o
mhd_mod.o     : extras_mod.o riemann_mod.o link_mod.o flux_mod.o divb_clean_mod.o \
                interface_mod.o
