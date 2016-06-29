# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
amr_commons.o           : hydro_parameters.o amr_parameters.o omp_timer_mod.o
amr_parameters          :
hydro_parameters.o      : amr_parameters.o
riemann_mod.o           : hydro_parameters.o amr_parameters.o
mhd_mod.o               : patch_mod.o riemann_mod.o hydro_parameters.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o
