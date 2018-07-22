# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/stagger2/Makefile


mhd_mod.o               : io_mod.o trace_mod.o scalar_mod.o vector_mod.o force_mod.o gpatch_mod.o scaling_mod.o units_mod.o download_mod.o \
                          $(STAGGER) eos_mod.o index_mod.o timestep_mod.o vector_ops.o kinds_mod.o non_ideal_mod.o viscosity_mod.o
force_mod.o             : io_mod.o trace_mod.o vector_mod.o mesh_mod.o index_mod.o vector_ops.o scalar_mod.o
viscosity_mod.o		: io_mod.o trace_mod.o mesh_mod.o force_mod.o
