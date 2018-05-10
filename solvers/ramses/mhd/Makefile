# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/ramses/Makefile

DEBUG += -DCHECK_NAN

# Dependencies
amr_parameters.o        :
hydro_parameters.o      : amr_parameters.o
amr_commons.o           : amr_parameters.o
hydro_commons.o         : amr_parameters.o hydro_parameters.o
units.o                 : amr_commons.o

riemann_mod.o           : hydro_parameters.o hydro_commons.o amr_parameters.o io_mod.o trace_mod.o
mhd_mod.o               : io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o patch_mod.o riemann_mod.o \
                          hydro_parameters.o link_mod.o download_mod.o gpatch_mod.o index_mod.o \
			  force_mod.o non_ideal_mod.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o hydro_parameters.o
index_mod.o             :
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o rt_mod.o
force_mod.o             : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o
non_ideal_mod.o         : io_mod.o trace_mod.o $(RAMSES) gpatch_mod.o link_mod.o force_mod.o
