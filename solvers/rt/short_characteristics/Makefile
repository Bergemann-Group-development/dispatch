VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)

rt_mod.o            : io_mod.o trace_mod.o radau_mod.o index_mod.o link_mod.o \
                      eos_mod.o gpatch_mod.o math_mod.o bits_mod.o rt_integral_mod.o \
                      scalar_mod.o timestep_mod.o rt_boundaries_mod.o aux_mod.o
rt_integral_mod.o   : io_mod.o trace_mod.o mesh_mod.o omp_mod.o
rt_solver_mod.o     : rt_mod.o
radau_mod.o	        : 
rt_boundaries_mod.o : trace_mod.o kinds_mod.o bits_mod.o mesh_mod.o eos_mod.o \
                      units_mod.o scaling_mod.o
rt_nbors_mod.o      : solver_mod.o
