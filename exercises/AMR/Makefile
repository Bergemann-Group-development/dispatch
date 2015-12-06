# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target
#PDE = advection
PDE = sedov
default: $(PDE).x

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


OBJ = \
  omp_timer_mod.o upwind_mod.o central_mod.o derivatives_mod.o interpolation_mod.o \
  io_mod.o patch_list_mod.o trace_mod.o omp_mod.o \
  scalar_mod.o vector_mod.o tensor_mod.o stagger_mod.o \
  bits_mod.o $(PDE)_mod.o patches_mod.o refine_mod.o timestep_mod.o

# Dependencies
patches_mod.o           : trace_mod.o io_mod.o
advection_mod.o         : patches_mod.o scalar_mod.o vector_mod.o tensor_mod.o trace_mod.o
patch_list_mod.o        : io_mod.o timestep_mod.o refine_mod.o
derivatives_mod.o       : upwind_mod.o
central_mod.o           : io_mod.o trace_mod.o
upwind_mod.o            : central_mod.o
timestep_mod.o          : io_mod.o $(PDE)_mod.o
refine_mod.o            : timestep_mod.o
patches_mod.o           : bits_mod.o central_mod.o
experiment.o            : io_mod.o patch_list_mod.o
interpolation_mod.o     : io_mod.o trace_mod.o patches_mod.o

stagger_mod.o           : io_mod.o trace_mod.o
sedov_mod.o             : scalar_mod.o vector_mod.o tensor_mod.o patches_mod.o stagger_mod.o

experiment.x: $(OBJ)
$(PDE).x: experiment.x
	cp -p experiment.x $(PDE).x
