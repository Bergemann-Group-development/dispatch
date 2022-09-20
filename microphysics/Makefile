# Set search path
VPATH := $(VPATH):$(TOP)/microphysics

# Source the EOS and opacity Makefiles
include $(TOP)/microphysics/eos/Makefile
include $(TOP)/data/opacities/Makefile

# Dependencies
itable_mod.o               : io_mod.o trace_mod.o kinds_mod.o
table_mod.o                : io_mod.o io_unit_mod.o scaling_mod.o trace_mod.o
microphysics_mod.o         : io_mod.o io_unit_mod.o trace_mod.o scaling_mod.o units_mod.o eos_mod.o opacity_mod.o
microphysics_factory_mod.o : io_mod.o io_unit_mod.o trace_mod.o scaling_mod.o units_mod.o eos_mod.o opacity_mod.o \
                             radiative_transfer_mod.o microphysics_mod.o
units_mod.o                : io_unit_mod.o trace_mod.o math_mod.o omp_mod.o io_mod.o
standard_scaling_mod.o     : io_mod.o units_mod.o shared_mod.o
scaling_mod.o              : standard_scaling_mod.o
bifrost_scaling_mod.o      : io_mod.o units_mod.o math_mod.o
math_mod.o                 :
