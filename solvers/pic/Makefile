ifeq ($(PREC),SINGLE)
  DEFS := $(DEFS) -DPREC=4
else
  DEFS := $(DEFS) -DPREC=8
endif


sinclude $(TOP)/solvers/pic/stagger/Makefile
sinclude $(TOP)/solvers/pic/Makefile.dep

VPATH := $(VPATH):$(TOP)/solvers/pic

pic_mod.o: io_mod.o mpi_mod.o trace_mod.o mhd_mod.o pic_stagger_$(STAGGER).o \
           grid.o species.o $(MPI).o maxwell_sources.o random.o \
           particleoperations.o textio.o

$(MPI).o:  stat.o textio.o

initial_pic_mod.o : io_mod.o trace_mod.o mesh_mod.o random_mod.o \
           initial_mod.o pic_mod.o

initial.o: io_mod.o trace_mod.o mesh_mod.o random_mod.o \
           initial_mod.o pic_mod.o index_mod.o
solver_mod.o: mhd_mod.o
