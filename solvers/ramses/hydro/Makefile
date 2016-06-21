# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
amr_commons.o           : hydro_parameters.o amr_parameters.o omp_timer_mod.o
godunov_fine.o          : hydro_commons.o amr_commons.o
godunov_utils.o         : hydro_parameters.o amr_parameters.o hydro_commons.o
hydro_commons.o         : amr_parameters.o hydro_parameters.o
hydro_parameters.o      : amr_parameters.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o
mhd_mod.o               : patch_mod.o godunov_fine.o godunov_utils.o umuscl.o \
                          units.o
umuscl.o                : hydro_commons.o amr_parameters.o hydro_parameters.o
units.o                 : amr_commons.o hydro_commons.o
