# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target
PDE = BE
EXEC= experiment.x

# Common configuration
TOP = ../..
sinclude $(TOP)/config/Makefile
TOP = ../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../../../..
sinclude $(TOP)/config/Makefile
VPATH:=$(VPATH):../

OBJ = \
  mpi_mod.o omp_mod.o omp_timer_mod.o io_mod.o trace_mod.o stagger_mod.o \
  bits_mod.o patches_mod.o timestep_mod.o patch_list_mod.o refine_mod.o \
  $(PDE)_mod.o hd_mod.o ghd_mod.o laplace_mod.o scalar_mod.o vector_mod.o tensor_mod.o

# Dependencies
$(PDE)_mod.o            : omp_timer_mod.o io_mod.o trace_mod.o patches_mod.o scalar_mod.o laplace_mod.o hd_mod.o ghd_mod.o
patches_mod.o           : trace_mod.o io_mod.o bits_mod.o
laplace_mod.o           : io_mod.o trace_mod.o patches_mod.o mpi_mod.o
hd_mod.o                : io_mod.o trace_mod.o scalar_mod.o vector_mod.o tensor_mod.o
ghd_mod.o               : hd_mod.o

# Dependencies
patch_list_mod.o        : io_mod.o timestep_mod.o refine_mod.o
central_mod.o           : io_mod.o trace_mod.o
timestep_mod.o          : io_mod.o $(PDE)_mod.o
refine_mod.o            : timestep_mod.o
patches_mod.o           : io_mod.o trace_mod.o bits_mod.o
experiment.o            : io_mod.o patch_list_mod.o
interpolation_mod.o     : io_mod.o trace_mod.o patches_mod.o
stagger_mod.o           : io_mod.o trace_mod.o
sedov_mod.o             : scalar_mod.o vector_mod.o tensor_mod.o patches_mod.o stagger_mod.o

$(EXEC): $(OBJ)
