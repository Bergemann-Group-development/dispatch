# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target
PDE = poisson
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

OBJ = \
  omp_mod.o omp_timer_mod.o io_mod.o trace_mod.o \
  $(PDE)_mod.o patches_mod.o scalar_mod.o laplace_mod.o

# Dependencies
$(PDE)_mod.o            : io_mod.o trace_mod.o patches_mod.o scalar_mod.o laplace_mod.o
patches_mod.o           : trace_mod.o io_mod.o bits_mod.o

$(EXEC): $(OBJ)
