VPATH := $(VPATH):$(TOP)/solvers/rt/$(RT_SOLVER)

rt_sc_mod.o	: io_mod.o trace_mod.o timer_mod.o link_mod.o patch_mod.o
rt_patch_mod.o  : io_mod.o trace_mod.o radau_mod.o solver_mod.o index_mod.o link_mod.o rt_sc_mod.o eos_mod.o
rt_mod.o        : io_mod.o trace_mod.o radau_mod.o index_mod.o link_mod.o eos_mod.o solver_mod.o \
                  math_mod.o
rt_solver_mod.o : rt_mod.o
rt_setup.o      : rt_mod.o task_list_mod.o data_io_mod.o
radau_mod.o	: 
