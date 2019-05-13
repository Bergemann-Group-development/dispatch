# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER):$(TOP)/solvers

include $(TOP)/solvers/ramses/Makefile

# Dependencies
amr_parameters          :
hydro_parameters.o      : amr_parameters.o
hydro_commons.o         : hydro_parameters.o
riemann_mod.o           : hydro_parameters.o hydro_commons.o amr_parameters.o io_mod.o trace_mod.o
hd_mod.o                : io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o extras_mod.o riemann_mod.o \
                          hydro_parameters.o link_mod.o download_mod.o eos_mod.o index_mod.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o hydro_parameters.o \
                          index_mod.o hydro_commons.o
force_mod.o             : io_unit_mod.o trace_mod.o io_mod.o mesh_mod.o random_mod.o scaling_mod.o
rhd_mod.o               : io_mod.o trace_mod.o hd_mod.o rt_mod.o
solver_mod.o            : io_mod.o trace_mod.o kinds_mod.o patch_mod.o hd_mod.o stagger_2nd.o
