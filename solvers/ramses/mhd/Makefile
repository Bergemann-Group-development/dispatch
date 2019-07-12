# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER):$(TOP)/solvers/$(SOLVER)/..

include $(TOP)/solvers/ramses/Makefile

# Dependencies
amr_parameters.o        :
hydro_parameters.o      : amr_parameters.o
amr_commons.o           : amr_parameters.o
hydro_commons.o         : amr_parameters.o hydro_parameters.o
units.o                 : amr_commons.o

riemann_mod.o           : hydro_parameters.o hydro_commons.o amr_parameters.o io_mod.o trace_mod.o
mhd_mod.o               : io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o patch_mod.o riemann_mod.o \
                          hydro_parameters.o link_mod.o download_mod.o extras_mod.o index_mod.o \
			  non_ideal_mod.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o hydro_parameters.o
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o rt_mod.o
force_mod.o             : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o
non_ideal_mod.o         : io_mod.o trace_mod.o link_mod.o force_mod.o
solver_mod.o            : io_mod.o trace_mod.o kinds_mod.o patch_mod.o mhd_mod.o validate_mod.o \
                          task_mod.o stagger_2nd.o
