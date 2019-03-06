VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)

rt_mod.o          : io_mod.o trace_mod.o radau_mod.o index_mod.o link_mod.o \
                    eos_mod.o solver_mod.o math_mod.o bits_mod.o rt_integral_mod.o \
                    scalar_mod.o timestep_mod.o
rt_integral_mod.o : io_mod.o trace_mod.o mesh_mod.o omp_mod.o
rt_solver_mod.o   : rt_mod.o
rt_setup.o        : rt_mod.o task_list_mod.o data_io_mod.o
radau_mod.o	      : 
