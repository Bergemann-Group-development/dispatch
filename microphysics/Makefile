# Set search path
VPATH := $(VPATH):$(TOP)/microphysics

sinclude $(TOP)/microphysics/eos/Makefile
sinclude $(TOP)/microphysics/opacity/Makefile

# Dependencies
microphysics_mod.o      : opacities_mod.o scaling_mod.o io_mod.o source_function_mod.o cgs_mod.o
$(EOS)_eos_mod.o        : io_mod.o io_unit_mod.o mpi_mod.o scaling_mod.o
source_function_mod.o   : io_mod.o trace_mod.o scaling_mod.o cgs_mod.o $(EOS)_eos_mod.o
