VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

# Dependencies
amr_commons.o:      amr_parameters.o omp_timer_mod.o
hydro_commons.o:    hydro_parameters.o
hydro_parameters.o: amr_parameters.o
poisson_commons.o:  poisson_parameters.o amr_commons.o amr_parameters.o
umuscl.o:           hydro_commons.o godunov_utils.o
godunov_fine.o:     poisson_commons.o hydro_commons.o utilities.o units.o
godunov_utils.o:    hydro_parameters.o amr_commons.o
courant_fine.o:     poisson_commons.o hydro_commons.o
uplmde.o:           hydro_commons.o
units.o:            hydro_commons.o amr_commons.o
utilities.o:        amr_commons.o

method_mod.o:       mesh_mod.o
utils_mod.o:        task_mod.o list_mod.o
interpol_mod.o:     task_mod.o

initial_mod.o:      io_mod.o trace_mod.o random_mod.o
mhd_mod.o:          patch_mod.o umuscl.o godunov_fine.o courant_fine.o
