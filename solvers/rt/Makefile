VPATH := $(VPATH):$(TOP)/solvers/rt

sinclude $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile

rt_setup_mod.o      : rt_mod.o task_list_mod.o data_io_mod.o io_unit_mod.o omp_mod.o \
                      link_mod.o task_mod.o patch_mod.o solver_mod.o trace_mod.o     \
                      bits_mod.o nbors_mod.o
