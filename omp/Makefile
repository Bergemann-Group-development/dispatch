LOCKS ?= omp

VPATH := $(VPATH):$(TOP)/$(LOCKS):$(TOP)/omp

omp_mod.o:	  io_unit_mod.o
omp_timer_mod.o:  omp_mod.o io_unit_mod.o
omp_lock_mod.o:   omp_mod.o omp_timer_mod.o shared_mod.o
read_lock_mod.o:  omp_lock_mod.o
