# Set search path
VPATH := $(VPATH):$(TOP)/microphysics

# Source the EOS and opacity Makefiles
include $(TOP)/microphysics/eos/Makefile
include $(TOP)/data/opacities/Makefile

# Dependencies
microphysics_mod.o         : io_mod.o io_unit_mod.o trace_mod.o scaling_mod.o units_mod.o eos_mod.o opacity_mod.o
microphysics_factory_mod.o : io_mod.o io_unit_mod.o trace_mod.o scaling_mod.o units_mod.o eos_mod.o opacity_mod.o \
                             radiative_transfer_mod.o microphysics_mod.o
