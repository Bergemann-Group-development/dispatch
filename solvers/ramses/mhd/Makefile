# Vpath
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)/..

include $(TOP)/solvers/ramses/Makefile

# Dependencies

riemann_mod.o           : parameters_mod.o io_mod.o trace_mod.o timer_mod.o eos_mod.o
mhd_mod.o               : parameters_mod.o io_unit_mod.o trace_mod.o mpi_mod.o omp_mod.o \
                          patch_mod.o riemann_mod.o link_mod.o download_mod.o extras_mod.o \
                          index_mod.o non_ideal_mod.o flux_mod.o timer_mod.o \
                          divb_clean_mod.o
initial_mod.o           : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o
rmhd_mod.o              : io_mod.o trace_mod.o mhd_mod.o rt_mod.o
force_mod.o             : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o
non_ideal_mod.o         : io_mod.o trace_mod.o link_mod.o force_mod.o
solver_mod.o            : io_mod.o trace_mod.o kinds_mod.o patch_mod.o mhd_mod.o validate_mod.o \
                          task_mod.o stagger_2nd.o
