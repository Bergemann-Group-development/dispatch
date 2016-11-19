# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
amr_parameters          :
hydro_parameters.o      : amr_parameters.o
riemann_mod.o           : hydro_parameters.o amr_parameters.o io_mod.o
mhd_mod.o               : io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o patch_mod.o riemann_mod.o hydro_parameters.o link_mod.o download_mod.o eos_mod.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o hydro_parameters.o
