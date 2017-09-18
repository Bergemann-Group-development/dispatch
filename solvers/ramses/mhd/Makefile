# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/ramses/Makefile

DEBUG += -DCHECK_NAN

# Dependencies
amr_parameters.o        :
pm_parameters.o         : amr_parameters.o
hydro_parameters.o      : amr_parameters.o
poisson_parameters.o    : amr_parameters.o
amr_commons.o           : amr_parameters.o
hydro_commons.o         : amr_parameters.o hydro_parameters.o
pm_commons.o            : amr_parameters.o pm_parameters.o random.o
poisson_commons.o       : amr_commons.o poisson_parameters.o
units.o                 : amr_commons.o
write_screen.o          : amr_commons.o hydro_commons.o poisson_commons.o pm_commons.o
random.o                : amr_commons.o

riemann_mod.o           : hydro_parameters.o hydro_commons.o amr_parameters.o io_mod.o trace_mod.o
mhd_mod.o               : io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o patch_mod.o riemann_mod.o \
                          hydro_parameters.o link_mod.o download_mod.o gpatch_mod.o index_mod.o \
			  output_amr.o output_hydro.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o hydro_parameters.o
index_mod.o             :
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o rt_mod.o
output_amr.o            : amr_commons.o hydro_commons.o pm_commons.o title.o write_screen.o units.o
output_hydro.o          : amr_commons.o hydro_commons.o title.o output_poisson.o #output_part.o
output_poisson.o        : amr_commons.o poisson_commons.o

#output_part.o           : amr_commons.o hydro_commons.o title.o output_poisson.o pm_commons.o \
#                          utilities.o particle_tree.o add_list.o remove_list.o nbors_utils.o
