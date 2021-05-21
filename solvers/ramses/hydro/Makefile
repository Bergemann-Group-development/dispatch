# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)/..

include $(TOP)/solvers/ramses/Makefile

# Dependencies
#amr_parameters          :
#hydro_parameters.o      : amr_parameters.o
#hydro_commons.o         : hydro_parameters.o
riemann_mod.o           : parameters_mod.o io_mod.o trace_mod.o eos_mod.o
hd_mod.o                : parameters_mod.o io_unit_mod.o trace_mod.o mpi_mod.o \
                          omp_mod.o extras_mod.o riemann_mod.o link_mod.o \
                          download_mod.o eos_mod.o index_mod.o
force_mod.o             : io_unit_mod.o trace_mod.o io_mod.o mesh_mod.o \
                          random_mod.o scaling_mod.o
rhd_mod.o               : io_mod.o trace_mod.o hd_mod.o rt_mod.o
solver_mod.o            : io_mod.o trace_mod.o kinds_mod.o patch_mod.o hd_mod.o \
                          stagger_2nd.o
