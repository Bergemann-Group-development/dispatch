ifeq ($(PREC),SINGLE)
  DEFS := $(DEFS) -DPREC=4
else
  DEFS := $(DEFS) -DPREC=8
endif

# To include the PPcode, PPCODE should point to the PPcode top directory
PPCODE ?= $(TOP)/solvers/pic/ppcode

BOUNDARY       ?= periodic
STAGGER_ORDER  ?= 2nd
INTERPOLATION  ?= tsc
PPMPI          ?= mpi

# Object file list (automatic generation is not useful here)
PICOBJ ?= params.o header.o debug.o grid.o units.o species.o dump.o stat.o textio.o \
  omp.o $(PPMPI).o distributions.o particleoperations.o io.o random.o fileoperations.o errors.o trace.o \
  pic_stagger_$(STAGGER_ORDER).o mesh_sources.o maxwell_sources.o math.o laplace.o energy_filter.o \
  maxwell_solver_$(STAGGER_ORDER).o load_balance.o sort.o init_ppcode.o bitoperations.o timer.o set_dt.o \
  move_particles.o celloperations.o interpolation_$(INTERPOLATION)_3d.o $(BOUNDARY).o \
  tests.o picmhd.o init_ppcode.o pic_mod.o

sinclude $(TOP)/solvers/pic/Makefile.dep

# Search paths for pic/PP module
VPATH:=\
 $(PPCODE):\
 $(PPCODE)/Boundaries:\
 $(PPCODE)/Fields:\
 $(PPCODE)/Init:\
 $(PPCODE)/Interpolation:\
 $(PPCODE)/Io:\
 $(PPCODE)/Mpi:\
 $(PPCODE)/Utilities:\
 $(PPCODE)/Particles:\
 $(TOP)/solvers/pic:\
 $(VPATH)

pic_particles_mod.o     : io_mod.o trace_mod.o patch_mod.o link_mod.o units.o dll_mod.o

pic_mod.o               : gpatch_mod.o link_mod.o io_mod.o trace_mod.o species.o eos_mod.o\
                          pic_stagger_$(STAGGER_ORDER).o particleoperations.o textio.o \
                          random.o maxwell_sources.o vector_ops.o kinds_mod.o vector_mod.o \
                          scalar_mod.o
pic_mod2.o              : gpatch_mod.o link_mod.o io_mod.o trace_mod.o species.o eos_mod.o\
                          pic_stagger_$(STAGGER_ORDER).o particleoperations.o textio.o \
                          random.o maxwell_sources.o vector_ops.o kinds_mod.o vector_mod.o scalar_mod.o

$(MPI).o:  stat.o textio.o

initial_pic_mod.o : io_mod.o trace_mod.o mesh_mod.o random_mod.o \
           initial_mod.o pic_mod.o

initial.o: io_mod.o trace_mod.o mesh_mod.o random_mod.o \
           initial_mod.o pic_mod.o index_mod.o
