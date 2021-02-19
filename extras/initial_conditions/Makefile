
INITIAL_CONDITIONS ?= tests
VPATH := $(TOP)/extras/initial_conditions/$(INITIAL_CONDITIONS):$(VPATH)
sinclude $(TOP)/extras/initial_conditions/$(INITIAL_CONDITIONS)/Makefile

# Generic dependencies
initial_conditions_mod.o : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o \
                           index_mod.o units_mod.o eos_mod.o gpatch_mod.o
